# Metabolite resampling methods

After MALDI-MSI coordinates are registered to the transcriptomic coordinate system, MSI feature intensities are resampled onto RNA spatial units.

## Primary method: area-overlap weighted resampling

Each RNA spatial unit and each MSI pixel are represented as square sampling regions in the registered coordinate system. For each RNA spatial unit, nearby MSI pixels are identified with a KD-tree. If one or more MSI pixels overlap the RNA spatial unit, their metabolite intensities are averaged using the geometric overlap area as weights.

This method is implemented in the original module-based script:

- `code/used_upsample_analysis_Module.py`
- `code/run_area_weighted_resampling.py`

A cleaned modular implementation is also provided in `code/jsrs/area_weighted_resampling.py` for code inspection and reuse.

## Alternative robustness methods

Two alternative methods are provided to assess whether downstream results depend strongly on the primary resampling strategy.

1. **Nearest-neighbour resampling**: each RNA spatial unit receives the metabolite profile of the closest registered MSI pixel.
2. **Inverse-distance weighted resampling**: each RNA spatial unit receives a weighted average of the `k` nearest MSI pixels. The default analysis uses `k=4` and `power=2`.

These methods are implemented in:

- `code/jsrs/alternative_resampling.py`
- `code/run_nearest_idw_resampling.py`

## Robustness analysis

The downstream robustness analysis compares cluster-level metabolite abundance, cluster similarity structure, log2FC concordance across cluster pairs, direction agreement and top-feature overlap across the three resampling methods.

Scripts are provided in `analysis/interpolation_robustness`.
