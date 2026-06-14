from __future__ import annotations

"""Initial rigid registration utilities for spatial point clouds.

The manuscript workflow used the transcriptomic coordinate system as the
reference and aligned metabolomic coordinates before JSRS resampling.  This file
provides a cleaned, self-contained version of the initial rigid-alignment stage:

1. Optional coordinate scaling/orientation adjustment is applied outside or
   before calling these functions.
2. Point clouds are rasterized into binary images for a coarse image-level
   rotation/translation search.
3. The image-based estimate is mapped back to point coordinates.
4. An optional point-cloud refinement minimizes nearest-neighbour distances.

The resulting coordinates can be used as an initialization for subsequent
SURF-based feature matching and image-based refinement.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from skimage.draw import disk
from skimage.transform import rotate


@dataclass
class RegistrationResult:
    """Result summary for the initial rigid registration stage."""

    image_angle_degrees: float
    image_translation_x: float
    image_translation_y: float
    image_loss: float
    image_success: bool
    point_refinement_angle_degrees: float
    point_refinement_translation_x: float
    point_refinement_translation_y: float
    point_refinement_loss: float
    point_refinement_success: bool
    transformed_points: np.ndarray

    def to_dict(self) -> dict:
        out = asdict(self)
        out.pop("transformed_points", None)
        return out


def scale_points_to_image(points: np.ndarray, image_size: int = 512) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scale coordinates into image-pixel space.

    Returns the scaled points together with the original minimum and span so the
    same coordinate system can be reconstructed later.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must be an N x 2 coordinate array.")
    min_xy = points.min(axis=0)
    span_xy = np.ptp(points, axis=0) + 1e-12
    image_points = (points - min_xy) / span_xy * (image_size - 1)
    return image_points, min_xy, span_xy


def unscale_points_from_image(image_points: np.ndarray, min_xy: np.ndarray, span_xy: np.ndarray, image_size: int = 512) -> np.ndarray:
    """Map image-pixel coordinates back to the original coordinate scale."""
    return np.asarray(image_points, dtype=float) / (image_size - 1) * span_xy + min_xy


def points_to_image(points: np.ndarray, image_size: int = 512, point_radius: int = 2) -> np.ndarray:
    """Rasterize 2D point coordinates into a binary image."""
    image_points, _, _ = scale_points_to_image(points, image_size=image_size)
    image = np.zeros((image_size, image_size), dtype=np.float32)
    pixel_points = np.round(image_points).astype(int)
    for x, y in pixel_points:
        rr, cc = disk((y, x), radius=point_radius, shape=image.shape)
        image[rr, cc] = 1.0
    return image


def transform_points(points: np.ndarray, angle_degrees: float, translation: tuple[float, float]) -> np.ndarray:
    """Apply a 2D rotation and translation to point coordinates."""
    points = np.asarray(points, dtype=float)
    theta = np.deg2rad(angle_degrees)
    rotation = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]],
        dtype=float,
    )
    transformed = points @ rotation.T
    transformed[:, 0] += translation[0]
    transformed[:, 1] += translation[1]
    return transformed


def transform_image(image: np.ndarray, angle_degrees: float, translation: tuple[float, float]) -> np.ndarray:
    """Rotate and translate a binary image."""
    rotated = rotate(image, angle_degrees, resize=False, mode="constant", cval=0, preserve_range=True)
    translated = np.roll(rotated, int(round(translation[0])), axis=1)
    translated = np.roll(translated, int(round(translation[1])), axis=0)
    return translated


def compute_image_loss(reference_image: np.ndarray, moving_image: np.ndarray) -> float:
    """Return the squared-error loss between two images after cropping to common size."""
    min_shape = (
        min(reference_image.shape[0], moving_image.shape[0]),
        min(reference_image.shape[1], moving_image.shape[1]),
    )
    ref = reference_image[: min_shape[0], : min_shape[1]]
    mov = moving_image[: min_shape[0], : min_shape[1]]
    return float(np.sum((ref - mov) ** 2))


def image_loss_objective(params: np.ndarray, reference_image: np.ndarray, moving_image: np.ndarray) -> float:
    """Objective function for image-level rigid transformation."""
    angle, tx, ty = params
    transformed = transform_image(moving_image, angle, (tx, ty))
    return compute_image_loss(reference_image, transformed)


def point_cloud_loss(params: np.ndarray, reference_points: np.ndarray, moving_points: np.ndarray) -> float:
    """Nearest-neighbour point-cloud loss after transforming the moving points."""
    angle, tx, ty = params
    transformed = transform_points(moving_points, angle, (tx, ty))
    tree = cKDTree(np.asarray(reference_points, dtype=float))
    distances, _ = tree.query(transformed, k=1)
    return float(np.sum(distances**2))


def align_by_minimum_coordinate(moving_points: np.ndarray, reference_points: np.ndarray) -> np.ndarray:
    """Translate moving points so their lower-left bounding-box corner matches the reference."""
    moving = np.asarray(moving_points, dtype=float)
    reference = np.asarray(reference_points, dtype=float)
    offset = reference.min(axis=0) - moving.min(axis=0)
    return moving + offset


def rotate_orientation(points: np.ndarray, mode: str = "none") -> np.ndarray:
    """Apply simple orientation transforms used before rigid optimization.

    Parameters
    ----------
    mode:
        One of ``none``, ``clockwise90``, ``counterclockwise90``, ``rotate180``
        or ``flip_x``.
    """
    points = np.asarray(points, dtype=float).copy()
    if mode == "none":
        return points
    if mode == "clockwise90":
        return np.column_stack([points[:, 1], -points[:, 0]])
    if mode == "counterclockwise90":
        return np.column_stack([-points[:, 1], points[:, 0]])
    if mode == "rotate180":
        return -points
    if mode == "flip_x":
        out = points.copy()
        out[:, 0] = out[:, 0].max() - out[:, 0]
        return out
    raise ValueError("mode must be one of: none, clockwise90, counterclockwise90, rotate180, flip_x")


def estimate_image_based_transform(
    moving_points: np.ndarray,
    reference_points: np.ndarray,
    image_size: int = 512,
    point_radius: int = 2,
    initial: tuple[float, float, float] = (0.0, 0.0, 0.0),
):
    """Estimate a coarse rigid transform using rasterized point-cloud images."""
    moving_image_points, _, _ = scale_points_to_image(moving_points, image_size=image_size)
    reference_image_points, reference_min, reference_span = scale_points_to_image(reference_points, image_size=image_size)

    moving_image = points_to_image(moving_points, image_size=image_size, point_radius=point_radius)
    reference_image = points_to_image(reference_points, image_size=image_size, point_radius=point_radius)

    fit = minimize(
        image_loss_objective,
        np.array(initial, dtype=float),
        args=(reference_image, moving_image),
        method="Powell",
    )
    angle, tx, ty = fit.x
    transformed_image_points = transform_points(moving_image_points, angle, (tx, ty))
    transformed_points = unscale_points_from_image(transformed_image_points, reference_min, reference_span, image_size=image_size)

    return {
        "angle_degrees": float(angle),
        "translation_x": float(tx),
        "translation_y": float(ty),
        "loss": float(fit.fun),
        "success": bool(fit.success),
        "transformed_points": transformed_points,
        "reference_image": reference_image,
        "moving_image": moving_image,
    }


def _subsample_points(points: np.ndarray, max_points: int | None, seed: int = 1234) -> np.ndarray:
    """Return a deterministic subsample for optimization if needed."""
    points = np.asarray(points, dtype=float)
    if max_points is None or max_points <= 0 or points.shape[0] <= max_points:
        return points
    rng = np.random.default_rng(seed)
    idx = rng.choice(points.shape[0], size=max_points, replace=False)
    return points[np.sort(idx), :]


def refine_with_point_cloud_loss(
    moving_points: np.ndarray,
    reference_points: np.ndarray,
    initial: tuple[float, float, float] = (0.0, 0.0, 0.0),
    max_points: int | None = 5000,
    random_seed: int = 1234,
):
    """Refine transformed coordinates by minimizing nearest-neighbour distances.

    The optimization can use deterministic subsampling to keep public example
    runs fast on large point clouds.  The fitted transformation is then applied
    to the full moving-point array.
    """
    moving_opt = _subsample_points(moving_points, max_points=max_points, seed=random_seed)
    reference_opt = _subsample_points(reference_points, max_points=max_points, seed=random_seed + 1)
    fit = minimize(
        point_cloud_loss,
        np.array(initial, dtype=float),
        args=(reference_opt, moving_opt),
        method="Powell",
    )
    angle, tx, ty = fit.x
    transformed_points = transform_points(moving_points, angle, (tx, ty))
    return {
        "angle_degrees": float(angle),
        "translation_x": float(tx),
        "translation_y": float(ty),
        "loss": float(fit.fun),
        "success": bool(fit.success),
        "transformed_points": transformed_points,
    }


def estimate_initial_rigid_transform(
    source_points: np.ndarray,
    target_points: np.ndarray,
    image_size: int = 512,
    point_radius: int = 2,
    orientation: str = "none",
    coordinate_scale: float = 1.0,
    run_point_refinement: bool = True,
    max_points_for_refinement: int | None = 5000,
) -> RegistrationResult:
    """Run the cleaned initial rigid-registration pipeline.

    Parameters
    ----------
    source_points:
        Moving coordinates, typically the MSI coordinates before alignment.
    target_points:
        Reference coordinates, typically transcriptomic bin coordinates.
    image_size:
        Raster image size for the image-level optimization.
    point_radius:
        Radius used when rasterizing points.
    orientation:
        Optional simple orientation transform applied to source points before
        optimization.
    coordinate_scale:
        Optional multiplier applied to source coordinates before orientation
        transform. This is useful when raw coordinate units differ between
        modalities.
    run_point_refinement:
        Whether to run a second nearest-neighbour point-cloud refinement.
    max_points_for_refinement:
        Optional deterministic subsample size used during point-cloud refinement.
        The final transformation is still applied to all source points.
    """
    source = np.asarray(source_points, dtype=float) * float(coordinate_scale)
    target = np.asarray(target_points, dtype=float)
    source = rotate_orientation(source, mode=orientation)

    image_fit = estimate_image_based_transform(
        moving_points=source,
        reference_points=target,
        image_size=image_size,
        point_radius=point_radius,
    )
    after_image = align_by_minimum_coordinate(image_fit["transformed_points"], target)

    if run_point_refinement:
        point_fit = refine_with_point_cloud_loss(after_image, target, max_points=max_points_for_refinement)
        final_points = point_fit["transformed_points"]
    else:
        point_fit = {
            "angle_degrees": 0.0,
            "translation_x": 0.0,
            "translation_y": 0.0,
            "loss": point_cloud_loss(np.array([0.0, 0.0, 0.0]), target, after_image),
            "success": True,
            "transformed_points": after_image,
        }
        final_points = after_image

    return RegistrationResult(
        image_angle_degrees=float(image_fit["angle_degrees"]),
        image_translation_x=float(image_fit["translation_x"]),
        image_translation_y=float(image_fit["translation_y"]),
        image_loss=float(image_fit["loss"]),
        image_success=bool(image_fit["success"]),
        point_refinement_angle_degrees=float(point_fit["angle_degrees"]),
        point_refinement_translation_x=float(point_fit["translation_x"]),
        point_refinement_translation_y=float(point_fit["translation_y"]),
        point_refinement_loss=float(point_fit["loss"]),
        point_refinement_success=bool(point_fit["success"]),
        transformed_points=final_points,
    )


def plot_point_cloud_overlay(
    reference_points: np.ndarray,
    moving_points: np.ndarray,
    output_path: Optional[str | Path] = None,
    title: str = "Point-cloud overlay",
) -> None:
    """Plot reference and transformed moving point clouds in one coordinate system."""
    reference = np.asarray(reference_points, dtype=float)
    moving = np.asarray(moving_points, dtype=float)
    plt.figure(figsize=(7, 7))
    plt.scatter(reference[:, 0], reference[:, 1], s=2, alpha=0.35, label="reference")
    plt.scatter(moving[:, 0], moving[:, 1], s=2, alpha=0.35, label="moving transformed")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(markerscale=4)
    plt.tight_layout()
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=300)
    plt.close()
