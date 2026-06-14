from __future__ import annotations

"""Area-overlap weighted MSI-to-transcriptome resampling.

The primary JSRS resampling step treats each transcriptomic bin and MSI pixel as
an axis-aligned square in the registered coordinate system. Candidate MSI pixels
near a transcriptomic bin are found with ``scipy.spatial.cKDTree``. When multiple
MSI pixels overlap the transcriptomic bin, their intensities are averaged using
weights proportional to the geometric overlap area.

This implementation keeps the detailed behavior of the original analysis code
while removing local paths and private project-specific identifiers.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


@dataclass
class ResamplingDiagnostics:
    """Summary of an area-weighted resampling run."""

    n_rna_bins: int
    n_msi_pixels: int
    n_features: int
    n_bins_with_candidates: int
    n_bins_with_positive_overlap: int
    n_bins_fallback_to_nearest: int
    n_bins_without_signal: int
    percent_bins_with_signal: float
    transcriptomics_side: float
    metabolomics_side: float
    search_radius: float

    def to_dict(self) -> dict:
        return asdict(self)


def compute_overlap_area(
    tx: float,
    ty: float,
    mx: float,
    my: float,
    transcriptomics_side: float,
    metabolomics_side: float,
) -> float:
    """Compute the overlap area between one RNA bin square and one MSI pixel square.

    Parameters
    ----------
    tx, ty:
        Centre coordinate of the transcriptomic bin.
    mx, my:
        Centre coordinate of the MSI pixel.
    transcriptomics_side:
        Side length of the transcriptomic square in the registered coordinate
        system.
    metabolomics_side:
        Side length of the MSI square in the registered coordinate system.
    """
    x_overlap = max(
        0.0,
        min(tx + transcriptomics_side / 2.0, mx + metabolomics_side / 2.0)
        - max(tx - transcriptomics_side / 2.0, mx - metabolomics_side / 2.0),
    )
    y_overlap = max(
        0.0,
        min(ty + transcriptomics_side / 2.0, my + metabolomics_side / 2.0)
        - max(ty - transcriptomics_side / 2.0, my - metabolomics_side / 2.0),
    )
    return float(x_overlap * y_overlap)


def _validate_inputs(
    rna_coords: pd.DataFrame,
    msi_coords: pd.DataFrame,
    msi_matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_rna = {"bin_id", "x", "y"}
    required_msi = {"spot_id", "x", "y"}
    if not required_rna.issubset(rna_coords.columns):
        raise ValueError(f"rna_coords must contain columns: {sorted(required_rna)}")
    if not required_msi.issubset(msi_coords.columns):
        raise ValueError(f"msi_coords must contain columns: {sorted(required_msi)}")
    if msi_matrix.empty:
        raise ValueError("msi_matrix is empty")

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
    return rna, msi, matrix


def _weighted_profile_for_one_bin(
    tx: float,
    ty: float,
    candidate_indices: list[int],
    msi_xy: np.ndarray,
    msi_values: np.ndarray,
    transcriptomics_side: float,
    metabolomics_side: float,
    fallback_to_nearest: bool = True,
) -> tuple[np.ndarray, bool, bool]:
    """Return the area-weighted MSI profile for one transcriptomic bin.

    Returns
    -------
    values:
        Resampled feature vector.
    has_positive_overlap:
        Whether at least one candidate pixel had positive square-overlap area.
    used_nearest_fallback:
        Whether the nearest candidate was used because no positive overlap area
        was observed.
    """
    if len(candidate_indices) == 0:
        return np.full(msi_values.shape[1], np.nan), False, False

    # The original analysis used direct assignment when only one MSI pixel was
    # found. This avoids unnecessary numerical edge cases and preserves the
    # behavior of the analysis script used for the manuscript.
    if len(candidate_indices) == 1:
        return msi_values[candidate_indices[0]], True, False

    overlap_areas = np.array(
        [
            compute_overlap_area(
                tx=tx,
                ty=ty,
                mx=msi_xy[j, 0],
                my=msi_xy[j, 1],
                transcriptomics_side=transcriptomics_side,
                metabolomics_side=metabolomics_side,
            )
            for j in candidate_indices
        ],
        dtype=float,
    )

    total_overlap = float(np.nansum(overlap_areas))
    if total_overlap <= 0:
        if not fallback_to_nearest:
            return np.full(msi_values.shape[1], np.nan), False, False
        distances = np.sqrt((msi_xy[candidate_indices, 0] - tx) ** 2 + (msi_xy[candidate_indices, 1] - ty) ** 2)
        nearest_idx = candidate_indices[int(np.argmin(distances))]
        return msi_values[nearest_idx], False, True

    weights = overlap_areas / total_overlap
    resampled = np.nansum(msi_values[candidate_indices, :] * weights[:, None], axis=0)
    return resampled, True, False


def resample_area_weighted(
    rna_coords: pd.DataFrame,
    msi_coords: pd.DataFrame,
    msi_matrix: pd.DataFrame,
    transcriptomics_side: float = 30.0,
    metabolomics_side: float = 40.0,
    search_radius: float = 28.0,
    fallback_to_nearest: bool = True,
    return_diagnostics: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, ResamplingDiagnostics]:
    """Resample registered MSI intensities onto transcriptomic spatial units.

    Parameters
    ----------
    rna_coords:
        DataFrame with columns ``bin_id``, ``x`` and ``y``.
    msi_coords:
        DataFrame with columns ``spot_id``, ``x`` and ``y``. Coordinates must
        already be transformed into the transcriptomic coordinate system.
    msi_matrix:
        MSI intensity matrix indexed by ``spot_id`` and with columns representing
        m/z or annotated metabolite features.
    transcriptomics_side:
        Side length of a transcriptomic bin in the registered coordinate system.
    metabolomics_side:
        Side length of an MSI pixel in the registered coordinate system.
    search_radius:
        Radius used to query candidate MSI pixels around each transcriptomic bin.
    fallback_to_nearest:
        If candidate pixels are found but square overlap is zero, assign the
        nearest candidate instead of returning missing values.
    return_diagnostics:
        If ``True``, also return a ``ResamplingDiagnostics`` object.
    """
    if transcriptomics_side <= 0 or metabolomics_side <= 0 or search_radius <= 0:
        raise ValueError("transcriptomics_side, metabolomics_side and search_radius must be positive.")

    rna, msi, matrix = _validate_inputs(rna_coords, msi_coords, msi_matrix)
    feature_names = list(matrix.columns)
    msi_values = matrix.to_numpy(dtype=float)
    msi_xy = msi[["x", "y"]].to_numpy(dtype=float)
    rna_xy = rna[["x", "y"]].to_numpy(dtype=float)

    tree = cKDTree(msi_xy)
    out_values = np.full((rna.shape[0], len(feature_names)), np.nan, dtype=float)

    n_with_candidates = 0
    n_with_positive_overlap = 0
    n_fallback = 0

    for i, (tx, ty) in enumerate(rna_xy):
        candidate_indices = tree.query_ball_point([tx, ty], r=search_radius)
        if candidate_indices:
            n_with_candidates += 1
        profile, has_overlap, used_fallback = _weighted_profile_for_one_bin(
            tx=float(tx),
            ty=float(ty),
            candidate_indices=list(candidate_indices),
            msi_xy=msi_xy,
            msi_values=msi_values,
            transcriptomics_side=transcriptomics_side,
            metabolomics_side=metabolomics_side,
            fallback_to_nearest=fallback_to_nearest,
        )
        out_values[i, :] = profile
        n_with_positive_overlap += int(has_overlap)
        n_fallback += int(used_fallback)

    result = rna[["bin_id", "x", "y"]].copy()
    feature_table = pd.DataFrame(out_values, columns=feature_names, index=result.index)
    result = pd.concat([result, feature_table], axis=1)

    has_signal = feature_table.notna().any(axis=1)
    diagnostics = ResamplingDiagnostics(
        n_rna_bins=int(rna.shape[0]),
        n_msi_pixels=int(msi.shape[0]),
        n_features=int(len(feature_names)),
        n_bins_with_candidates=int(n_with_candidates),
        n_bins_with_positive_overlap=int(n_with_positive_overlap),
        n_bins_fallback_to_nearest=int(n_fallback),
        n_bins_without_signal=int((~has_signal).sum()),
        percent_bins_with_signal=float(has_signal.mean() * 100.0),
        transcriptomics_side=float(transcriptomics_side),
        metabolomics_side=float(metabolomics_side),
        search_radius=float(search_radius),
    )

    if return_diagnostics:
        return result, diagnostics
    return result


def batch_resample_area_weighted(
    sections: Iterable[str],
    rna_coord_template: str,
    msi_coord_template: str,
    msi_matrix_template: str,
    output_template: str,
    transcriptomics_side: float = 30.0,
    metabolomics_side: float = 40.0,
    search_radius: float = 28.0,
) -> list[dict]:
    """Run area-weighted resampling for multiple sections.

    Templates should contain ``{section}``, for example
    ``example_data/rna_coords_{section}.csv``.  This helper is mainly intended
    to document the multi-section workflow; the command-line example scripts use
    the same core function for one section at a time.
    """
    from .io_utils import read_msi_matrix, read_table, standardize_msi_coords, standardize_rna_coords, write_resampled_table

    summaries: list[dict] = []
    for section in sections:
        rna_coords = standardize_rna_coords(read_table(rna_coord_template.format(section=section)))
        msi_coords = standardize_msi_coords(read_table(msi_coord_template.format(section=section)))
        matrix = read_msi_matrix(msi_matrix_template.format(section=section), msi_coords["spot_id"])
        result, diagnostics = resample_area_weighted(
            rna_coords=rna_coords,
            msi_coords=msi_coords,
            msi_matrix=matrix,
            transcriptomics_side=transcriptomics_side,
            metabolomics_side=metabolomics_side,
            search_radius=search_radius,
            return_diagnostics=True,
        )
        write_resampled_table(result, output_template.format(section=section))
        d = diagnostics.to_dict()
        d["section"] = section
        summaries.append(d)
    return summaries
