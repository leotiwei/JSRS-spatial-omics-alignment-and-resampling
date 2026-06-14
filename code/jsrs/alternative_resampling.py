from __future__ import annotations

"""Alternative MSI-to-transcriptome resampling methods.

These methods were used for robustness analysis.  They operate on the same
registered coordinates and intensity matrices as the primary area-overlap
weighted JSRS resampling method.
"""

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


@dataclass
class AlternativeResamplingDiagnostics:
    method: str
    n_rna_bins: int
    n_msi_pixels: int
    n_features: int
    k: int | None
    power: float | None
    max_distance: float | None
    n_bins_with_signal: int
    percent_bins_with_signal: float

    def to_dict(self) -> dict:
        return asdict(self)


def _prepare_resampling_inputs(
    rna_coords: pd.DataFrame,
    msi_coords: pd.DataFrame,
    msi_matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[str]]:
    required_rna = {"bin_id", "x", "y"}
    required_msi = {"spot_id", "x", "y"}
    if not required_rna.issubset(rna_coords.columns):
        raise ValueError(f"rna_coords must contain columns: {sorted(required_rna)}")
    if not required_msi.issubset(msi_coords.columns):
        raise ValueError(f"msi_coords must contain columns: {sorted(required_msi)}")

    rna = rna_coords.copy()
    msi = msi_coords.copy()
    matrix = msi_matrix.copy()
    rna["bin_id"] = rna["bin_id"].astype(str)
    msi["spot_id"] = msi["spot_id"].astype(str)
    matrix.index = matrix.index.astype(str)

    common_spot_ids = [spot_id for spot_id in msi["spot_id"] if spot_id in matrix.index]
    if len(common_spot_ids) == 0:
        raise ValueError("No shared spot IDs between MSI coordinates and MSI matrix.")
    msi = msi.set_index("spot_id").loc[common_spot_ids].reset_index()
    matrix = matrix.loc[common_spot_ids]
    return rna, msi, matrix.to_numpy(dtype=float), list(matrix.columns)


def _assemble_output(rna_coords: pd.DataFrame, values: np.ndarray, feature_names: list[str]) -> pd.DataFrame:
    base = rna_coords[["bin_id", "x", "y"]].copy()
    feature_table = pd.DataFrame(values, columns=feature_names, index=base.index)
    return pd.concat([base, feature_table], axis=1)


def _diagnostics(
    method: str,
    values: np.ndarray,
    n_rna_bins: int,
    n_msi_pixels: int,
    n_features: int,
    k: int | None = None,
    power: float | None = None,
    max_distance: float | None = None,
) -> AlternativeResamplingDiagnostics:
    has_signal = ~np.isnan(values).all(axis=1)
    return AlternativeResamplingDiagnostics(
        method=method,
        n_rna_bins=int(n_rna_bins),
        n_msi_pixels=int(n_msi_pixels),
        n_features=int(n_features),
        k=k,
        power=power,
        max_distance=max_distance,
        n_bins_with_signal=int(has_signal.sum()),
        percent_bins_with_signal=float(has_signal.mean() * 100.0),
    )


def resample_nearest(
    rna_coords: pd.DataFrame,
    msi_coords: pd.DataFrame,
    msi_matrix: pd.DataFrame,
    max_distance: float | None = None,
    return_diagnostics: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, AlternativeResamplingDiagnostics]:
    """Assign each RNA bin the MSI profile of the nearest registered MSI pixel.

    Parameters
    ----------
    max_distance:
        Optional maximum accepted distance between an RNA bin and the nearest MSI
        pixel. If supplied, RNA bins farther than this distance are filled with
        missing values.
    """
    rna, msi, msi_values, feature_names = _prepare_resampling_inputs(rna_coords, msi_coords, msi_matrix)
    tree = cKDTree(msi[["x", "y"]].to_numpy(dtype=float))
    distances, indices = tree.query(rna[["x", "y"]].to_numpy(dtype=float), k=1)
    values = msi_values[indices, :].copy()

    if max_distance is not None:
        values[distances > max_distance, :] = np.nan

    result = _assemble_output(rna, values, feature_names)
    diag = _diagnostics(
        method="nearest",
        values=values,
        n_rna_bins=rna.shape[0],
        n_msi_pixels=msi.shape[0],
        n_features=len(feature_names),
        max_distance=max_distance,
    )
    if return_diagnostics:
        return result, diag
    return result


def resample_idw(
    rna_coords: pd.DataFrame,
    msi_coords: pd.DataFrame,
    msi_matrix: pd.DataFrame,
    k: int = 4,
    power: float = 2.0,
    eps: float = 1e-12,
    max_distance: float | None = None,
    return_diagnostics: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, AlternativeResamplingDiagnostics]:
    """Resample MSI intensities using inverse-distance weighting.

    For each RNA bin, the ``k`` nearest MSI pixels are found in the registered
    coordinate system. Feature intensities are averaged with weights
    ``1 / distance**power``.  Zero-distance matches are handled with ``eps`` so
    that co-located pixels receive dominant weight without causing division by
    zero.
    """
    if k <= 0:
        raise ValueError("k must be a positive integer.")
    if power <= 0:
        raise ValueError("power must be positive.")

    rna, msi, msi_values, feature_names = _prepare_resampling_inputs(rna_coords, msi_coords, msi_matrix)
    k_eff = min(int(k), msi.shape[0])
    tree = cKDTree(msi[["x", "y"]].to_numpy(dtype=float))
    distances, indices = tree.query(rna[["x", "y"]].to_numpy(dtype=float), k=k_eff)

    if k_eff == 1:
        distances = distances[:, None]
        indices = indices[:, None]

    weights = 1.0 / np.maximum(distances, eps) ** power
    if max_distance is not None:
        weights = np.where(distances <= max_distance, weights, 0.0)
    weight_sums = weights.sum(axis=1, keepdims=True)

    values = np.full((rna.shape[0], len(feature_names)), np.nan, dtype=float)
    valid_rows = weight_sums[:, 0] > 0
    normalized_weights = np.zeros_like(weights)
    normalized_weights[valid_rows, :] = weights[valid_rows, :] / weight_sums[valid_rows, :]
    values[valid_rows, :] = np.einsum("ij,ijk->ik", normalized_weights[valid_rows, :], msi_values[indices[valid_rows, :], :])

    result = _assemble_output(rna, values, feature_names)
    diag = _diagnostics(
        method="IDW",
        values=values,
        n_rna_bins=rna.shape[0],
        n_msi_pixels=msi.shape[0],
        n_features=len(feature_names),
        k=k_eff,
        power=float(power),
        max_distance=max_distance,
    )
    if return_diagnostics:
        return result, diag
    return result


def batch_resample_nearest_idw(
    sections: Iterable[str],
    rna_coord_template: str,
    msi_coord_template: str,
    msi_matrix_template: str,
    nearest_output_template: str,
    idw_output_template: str,
    k: int = 4,
    power: float = 2.0,
) -> list[dict]:
    """Run nearest-neighbour and IDW resampling for multiple sections."""
    from .io_utils import read_msi_matrix, read_table, standardize_msi_coords, standardize_rna_coords, write_resampled_table

    summaries: list[dict] = []
    for section in sections:
        rna_coords = standardize_rna_coords(read_table(rna_coord_template.format(section=section)))
        msi_coords = standardize_msi_coords(read_table(msi_coord_template.format(section=section)))
        matrix = read_msi_matrix(msi_matrix_template.format(section=section), msi_coords["spot_id"])

        nearest, nearest_diag = resample_nearest(rna_coords, msi_coords, matrix, return_diagnostics=True)
        idw, idw_diag = resample_idw(rna_coords, msi_coords, matrix, k=k, power=power, return_diagnostics=True)
        write_resampled_table(nearest, nearest_output_template.format(section=section))
        write_resampled_table(idw, idw_output_template.format(section=section))

        for diag in (nearest_diag, idw_diag):
            d = diag.to_dict()
            d["section"] = section
            summaries.append(d)
    return summaries
