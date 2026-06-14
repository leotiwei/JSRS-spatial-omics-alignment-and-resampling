# Software versions and external tools

The following tools were used in the analysis workflow. Exact versions may vary by analysis module and should be checked with the corresponding `sessionInfo()` or environment file when rerunning the scripts.

## Spatial transcriptomics

- Stereo-seq FF v1.3
- SAW 8.0.2
- Reference genome: GRCh38.105
- STAR aligner
- Seurat 5.1.0
- SeuratDisk 0.0.0.9021
- Scanpy 1.8.2

## Spatial metabolomics

- pyimzML 1.5.3
- HMDB accurate-mass annotation
- Scanpy 1.8.2

## Spatial registration and resampling

- Python 3.x
- NumPy
- pandas
- SciPy 1.13.1
- scikit-image
- SURF-based feature matching for registration refinement

## Downstream analyses

- Monocle 2.32.0
- CellChat
- SCENIC / AUCell
- R packages: Seurat, Matrix, data.table, dplyr, ggplot2, openxlsx, pheatmap
