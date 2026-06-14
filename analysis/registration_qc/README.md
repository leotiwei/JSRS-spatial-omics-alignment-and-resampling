# Registration QC

This directory contains scripts for evaluating RNA-to-MSI spatial registration accuracy using two complementary readouts:

1. Predefined anatomical landmark error before and after registration.
2. Tissue-mask overlap before and after registration.

The scripts start from processed coordinate tables and landmark tables. The example landmark table shows the required input format.

## Expected inputs

- `rna_coords_<section>.csv`: RNA / Stereo-seq coordinate table with columns `x` and `y`.
- `msi_coords_before_<section>.csv`: MSI coordinates before registration.
- `msi_coords_after_<section>.csv`: MSI coordinates after registration.
- `example_landmarks.tsv`: anatomical landmarks with RNA, before-registration MSI and after-registration MSI coordinates.

## Main output

- Landmark error before and after registration.
- Dice and Jaccard mask-overlap summaries.
- Overlay plots for visual inspection.
