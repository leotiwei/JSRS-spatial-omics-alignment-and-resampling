#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from jsrs.initial_rigid_registration import estimate_initial_rigid_transform, plot_point_cloud_overlay
from jsrs.io_utils import read_table, require_coordinate_columns, write_summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Estimate an initial rigid transform between two spatial point clouds.")
    parser.add_argument("--source-coords", required=True, help="Moving coordinate table, typically MSI before alignment.")
    parser.add_argument("--target-coords", required=True, help="Reference coordinate table, typically RNA coordinates.")
    parser.add_argument("--output", required=True, help="Output CSV containing transformed source coordinates.")
    parser.add_argument("--summary-output", default=None, help="Optional TSV file for registration parameters.")
    parser.add_argument("--plot-output", default=None, help="Optional PNG/PDF overlay plot path.")
    parser.add_argument("--image-size", type=int, default=256, help="Raster image size for image-based optimization.")
    parser.add_argument("--point-radius", type=int, default=2, help="Point radius for rasterization.")
    parser.add_argument(
        "--orientation",
        default="none",
        choices=["none", "clockwise90", "counterclockwise90", "rotate180", "flip_x"],
        help="Optional orientation transform applied before optimization.",
    )
    parser.add_argument("--coordinate-scale", type=float, default=1.0, help="Optional multiplier for source coordinates.")
    parser.add_argument("--skip-point-refinement", action="store_true", help="Skip point-cloud refinement after image optimization.")
    parser.add_argument("--max-points-for-refinement", type=int, default=5000, help="Maximum points used for point-cloud refinement.")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    source = read_table(args.source_coords)
    target = read_table(args.target_coords)
    sx, sy = require_coordinate_columns(source, "source coordinate table")
    tx, ty = require_coordinate_columns(target, "target coordinate table")

    fit = estimate_initial_rigid_transform(
        source_points=source[[sx, sy]].to_numpy(),
        target_points=target[[tx, ty]].to_numpy(),
        image_size=args.image_size,
        point_radius=args.point_radius,
        orientation=args.orientation,
        coordinate_scale=args.coordinate_scale,
        run_point_refinement=not args.skip_point_refinement,
        max_points_for_refinement=args.max_points_for_refinement,
    )

    output = source.copy()
    output["x"] = fit.transformed_points[:, 0]
    output["y"] = fit.transformed_points[:, 1]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    summary = fit.to_dict()
    if args.summary_output is not None:
        write_summary(summary, args.summary_output)
    if args.plot_output is not None:
        plot_point_cloud_overlay(target[[tx, ty]].to_numpy(), fit.transformed_points, args.plot_output)

    print(json.dumps(summary, indent=2))
    print(f"Wrote transformed coordinates to: {args.output}")


if __name__ == "__main__":
    main()
