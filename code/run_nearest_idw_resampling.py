#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from jsrs.alternative_resampling import resample_idw, resample_nearest
from jsrs.io_utils import (
    read_msi_matrix,
    read_table,
    standardize_msi_coords,
    standardize_rna_coords,
    write_resampled_table,
    write_summary,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run nearest-neighbour and IDW MSI-to-RNA resampling.")
    parser.add_argument("--rna-coords", required=True, help="RNA coordinate CSV/TSV with x, y and optional bin ID.")
    parser.add_argument("--msi-coords", required=True, help="Registered MSI coordinate CSV/TSV with x, y and optional spot ID.")
    parser.add_argument("--msi-matrix", required=True, help="MSI spot-by-feature intensity matrix.")
    parser.add_argument("--output-nearest", required=True, help="Output CSV for nearest-neighbour resampling.")
    parser.add_argument("--output-idw", required=True, help="Output CSV for inverse-distance weighted resampling.")
    parser.add_argument("--summary-output", default=None, help="Optional TSV file for run diagnostics.")
    parser.add_argument("--idw-k", type=int, default=4, help="Number of neighbours for IDW.")
    parser.add_argument("--idw-power", type=float, default=2.0, help="Inverse-distance power for IDW.")
    parser.add_argument("--max-distance", type=float, default=None, help="Optional maximum accepted neighbour distance.")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    rna_coords = standardize_rna_coords(read_table(args.rna_coords))
    msi_coords = standardize_msi_coords(read_table(args.msi_coords))
    msi_matrix = read_msi_matrix(args.msi_matrix, msi_coords["spot_id"])

    nearest, nearest_diag = resample_nearest(
        rna_coords=rna_coords,
        msi_coords=msi_coords,
        msi_matrix=msi_matrix,
        max_distance=args.max_distance,
        return_diagnostics=True,
    )
    idw, idw_diag = resample_idw(
        rna_coords=rna_coords,
        msi_coords=msi_coords,
        msi_matrix=msi_matrix,
        k=args.idw_k,
        power=args.idw_power,
        max_distance=args.max_distance,
        return_diagnostics=True,
    )

    write_resampled_table(nearest, args.output_nearest)
    write_resampled_table(idw, args.output_idw)

    summary = {"nearest": nearest_diag.to_dict(), "idw": idw_diag.to_dict()}
    if args.summary_output is not None:
        flat = {f"nearest.{k}": v for k, v in nearest_diag.to_dict().items()}
        flat.update({f"idw.{k}": v for k, v in idw_diag.to_dict().items()})
        write_summary(flat, args.summary_output)

    print(json.dumps(summary, indent=2))
    print(f"Wrote nearest matrix to: {args.output_nearest}")
    print(f"Wrote IDW matrix to: {args.output_idw}")


if __name__ == "__main__":
    main()
