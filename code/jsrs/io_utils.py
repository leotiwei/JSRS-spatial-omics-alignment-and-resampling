from __future__ import annotations

"""Input/output utilities for JSRS example workflows.

The functions in this file intentionally accept several common column-name
conventions because exported coordinate and matrix files from spatial omics
software are not always identical.  The goal is to make the public scripts
inspectable and reusable without binding them to local file paths or private
intermediate object names.
"""

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


BIN_ID_CANDIDATES = (
    "bin_id", "barcode", "Barcode", "cell", "cell_id", "spot", "spot_id",
    "_index", "index", "Unnamed: 0",
)
MSI_ID_CANDIDATES = (
    "spot_id", "Spot_index", "pixel_id", "msi_id", "barcode", "Barcode",
    "cell", "cell_id", "_index", "index", "Unnamed: 0",
)
COORD_X_CANDIDATES = ("x", "X", "coord_x", "image_x", "aligned_x")
COORD_Y_CANDIDATES = ("y", "Y", "coord_y", "image_y", "aligned_y")


def read_table(path: str | Path, sep: Optional[str] = None) -> pd.DataFrame:
    """Read a CSV/TSV text table.

    Parameters
    ----------
    path:
        Input table path.
    sep:
        Optional delimiter. If not supplied, the delimiter is inferred from the
        file extension: TSV/TXT files are read as tab-separated, and other files
        are read as comma-separated.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if sep is None:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def first_present(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    """Return the first candidate column name present in a table."""
    column_list = list(columns)
    for candidate in candidates:
        if candidate in column_list:
            return candidate
    return None


def require_coordinate_columns(df: pd.DataFrame, table_name: str = "coordinate table") -> tuple[str, str]:
    """Find x/y coordinate columns and raise a readable error if absent."""
    x_col = first_present(df.columns, COORD_X_CANDIDATES)
    y_col = first_present(df.columns, COORD_Y_CANDIDATES)
    if x_col is None or y_col is None:
        raise ValueError(
            f"{table_name} must contain coordinate columns. Expected one of "
            f"{COORD_X_CANDIDATES} for x and one of {COORD_Y_CANDIDATES} for y."
        )
    return x_col, y_col


def standardize_coordinates(
    df: pd.DataFrame,
    id_candidates: Iterable[str],
    output_id_col: str,
    table_name: str,
    generated_prefix: str,
) -> pd.DataFrame:
    """Return a coordinate table with standardized ID, x and y columns."""
    x_col, y_col = require_coordinate_columns(df, table_name)
    id_col = first_present(df.columns, id_candidates)

    out = df.copy()
    if id_col is None:
        out[output_id_col] = [f"{generated_prefix}_{i}" for i in range(out.shape[0])]
    else:
        out[output_id_col] = out[id_col].astype(str)

    out["x"] = pd.to_numeric(out[x_col], errors="raise")
    out["y"] = pd.to_numeric(out[y_col], errors="raise")
    return out[[output_id_col, "x", "y"]]


def standardize_rna_coords(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize transcriptomic spatial coordinates.

    Returns
    -------
    pandas.DataFrame
        Columns: ``bin_id``, ``x`` and ``y``.
    """
    return standardize_coordinates(
        df=df,
        id_candidates=BIN_ID_CANDIDATES,
        output_id_col="bin_id",
        table_name="RNA coordinate table",
        generated_prefix="bin",
    )


def standardize_msi_coords(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize registered MSI spatial coordinates.

    Returns
    -------
    pandas.DataFrame
        Columns: ``spot_id``, ``x`` and ``y``.
    """
    return standardize_coordinates(
        df=df,
        id_candidates=MSI_ID_CANDIDATES,
        output_id_col="spot_id",
        table_name="MSI coordinate table",
        generated_prefix="spot",
    )


def read_msi_matrix(path: str | Path, spot_ids: Iterable[str]) -> pd.DataFrame:
    """Read an MSI intensity matrix and return a spot-by-feature table.

    The function handles three common formats:

    1. The first column contains MSI spot identifiers and the remaining columns
       are m/z or metabolite features.
    2. The table is transposed, with features in rows and spot identifiers as
       columns. In this case the matrix is transposed automatically.
    3. No spot identifier is present, but the row count matches the coordinate
       table. In this fallback case, row order is assumed to match coordinates.

    All feature columns are coerced to numeric values; non-numeric values become
    ``NaN`` so that downstream averaging can skip or propagate missing values.
    """
    matrix = read_table(path)
    spot_ids = [str(x) for x in spot_ids]
    spot_set = set(spot_ids)

    if matrix.empty:
        raise ValueError(f"MSI matrix is empty: {path}")

    first_col = matrix.columns[0]
    first_values = matrix[first_col].astype(str).tolist()
    row_id_overlap = len(set(first_values).intersection(spot_set))
    col_id_overlap = len(set(map(str, matrix.columns)).intersection(spot_set))

    if row_id_overlap > 0:
        matrix = matrix.set_index(first_col)
        matrix.index = matrix.index.astype(str)
    elif col_id_overlap > 0:
        feature_col = first_col
        matrix = matrix.set_index(feature_col).T
        matrix.index = matrix.index.astype(str)
    else:
        if matrix.shape[0] != len(spot_ids):
            raise ValueError(
                "Could not match MSI matrix rows or columns to MSI coordinates. "
                "Provide spot identifiers in the first column or as column names, "
                "or provide a matrix whose row order exactly matches the MSI coordinate table."
            )
        matrix.index = spot_ids

    available = [spot_id for spot_id in spot_ids if spot_id in matrix.index]
    if len(available) == 0:
        raise ValueError("No MSI spot IDs from the coordinate table were found in the matrix.")
    matrix = matrix.loc[available]

    # Drop purely coordinate-like columns if they were accidentally included in
    # the intensity matrix. Keeping them would incorrectly treat coordinates as
    # metabolite features.
    drop_cols = [c for c in matrix.columns if str(c) in {"x", "y", "X", "Y", "spot_id", "Spot_index"}]
    if drop_cols:
        matrix = matrix.drop(columns=drop_cols)

    for column in matrix.columns:
        matrix[column] = pd.to_numeric(matrix[column], errors="coerce")

    return matrix


def summarize_matrix(matrix: pd.DataFrame) -> dict[str, float | int]:
    """Return basic feature-matrix QC statistics."""
    values = matrix.to_numpy(dtype=float)
    detected_per_row = np.sum(~np.isnan(values), axis=1)
    return {
        "n_rows": int(matrix.shape[0]),
        "n_features": int(matrix.shape[1]),
        "median_detected_features_per_row": float(np.median(detected_per_row)) if matrix.shape[0] else 0.0,
        "missing_fraction": float(np.mean(np.isnan(values))) if values.size else 0.0,
    }


def write_resampled_table(result: pd.DataFrame, output_path: str | Path) -> None:
    """Write a resampled feature table and create parent folders if needed."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)


def write_summary(summary: dict, output_path: str | Path) -> None:
    """Write a dictionary summary as a two-column TSV file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"metric": list(summary.keys()), "value": list(summary.values())}).to_csv(
        output_path, sep="\t", index=False
    )
