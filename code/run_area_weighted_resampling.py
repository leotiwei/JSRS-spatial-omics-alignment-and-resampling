#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run area-overlap weighted metabolite resampling.

This wrapper follows the original module-based usage pattern used for the
area-weighted resampling step. The core implementation is kept in
`used_upsample_analysis_Module.py`.

The example uses the small demonstration files in `example_data/` and writes
outputs to `outputs/area_weighted_resampling/`.
"""

import os
import sys
from pathlib import Path

# Add the current code directory so the original module can be imported.
CURRENT_DIR = Path(__file__).resolve().parent
REPO_DIR = CURRENT_DIR.parent
sys.path.append(str(CURRENT_DIR))

from used_upsample_analysis_Module import upsample_and_analyze


def main() -> None:
    example_dir = REPO_DIR / "example_data"
    output_dir = REPO_DIR / "outputs" / "area_weighted_resampling"
    output_dir.mkdir(parents=True, exist_ok=True)

    upsample_and_analyze(
        P_df_path=str(example_dir / "rna_coords_cs12_3.csv"),
        Q_df_path=str(example_dir / "msi_coords_aligned_cs12_3.csv"),
        meta_df_path=str(example_dir / "metabolite_matrix_cs12_3.csv"),
        transcriptomics_side=30,
        metabolomics_side=40,
        search_radius=28,
        output_dir=str(output_dir),
        output_csv="cs12_3_combined_multiomics_data_radius28.csv",
        log_file="cs12_3_process_log.txt",
    )


if __name__ == "__main__":
    os.environ.setdefault("MPLBACKEND", "Agg")
    main()
