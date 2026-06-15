# Code folder

This folder contains the Python scripts and modules used for JSRS-related coordinate alignment and metabolite resampling. The scripts are provided to document the study-specific custom steps, while mature algorithms and established software packages are referenced rather than reimplemented.

## Folder contents

```text
code/
  README.md
  run_initial_rigid_registration.py
  run_area_weighted_resampling.py
  run_nearest_idw_resampling.py
  used_upsample_analysis_Module.py
  jsrs/
    __init__.py
    io_utils.py
    initial_rigid_registration.py
    area_weighted_resampling.py
    alternative_resampling.py
```

### `run_initial_rigid_registration.py`

Example script for the initial coordinate-alignment step. This script calls functions in `jsrs/initial_rigid_registration.py` to estimate an approximate rigid transformation between transcriptomic spatial coordinates and metabolomic spatial coordinates.

This initial alignment code was developed for the study and is included in this repository. It uses coordinate-derived raster images and point-cloud representations to estimate approximate rotation, translation and scaling parameters. The optimization step uses `scipy.optimize.minimize` with the Powell method.

### `jsrs/initial_rigid_registration.py`

Reusable functions for the initial alignment step, including:

- coordinate loading and normalization;
- coordinate scaling and orientation transformation;
- rasterization of spatial point clouds;
- image-level objective calculation;
- rigid transformation optimization;
- point-cloud refinement;
- diagnostic overlay plotting.

This module documents the custom initial alignment procedure used before image-based refinement.

### Further image-based refinement

After the initial coordinate transformation, further fine alignment was performed using a standard SURF-based feature-matching workflow on section images or coordinate-derived reference images. SURF (Speeded-Up Robust Features) is an established local feature detector and descriptor for image matching and registration. In this study, the custom initial alignment code was used to bring the two coordinate systems into approximate register, and the subsequent SURF-based feature matching provided image-level refinement before downstream metabolite resampling.

The SURF algorithm itself was not reimplemented in this repository. Instead, the workflow used the established SURF feature-detection and matching strategy, as described in the manuscript Methods and in the original SURF publications. OpenCV provides SURF functionality through `cv.xfeatures2d.SURF_create` in the OpenCV contrib modules.

### `run_area_weighted_resampling.py`

Wrapper script for the primary area-overlap weighted metabolite resampling step. This script follows the original project structure and calls `used_upsample_analysis_Module.py`.

### `used_upsample_analysis_Module.py`

Original module used for area-overlap weighted resampling of metabolite intensities from registered MSI pixels onto transcriptomic spatial bins. This module is retained to document the primary JSRS resampling implementation used in the analysis.

### `run_nearest_idw_resampling.py`

Example script for alternative metabolite resampling strategies used in robustness analysis. It calls functions in `jsrs/alternative_resampling.py`.

### `jsrs/alternative_resampling.py`

Reusable functions for nearest-neighbour and inverse-distance weighted resampling. These alternative strategies were used to test whether downstream metabolomic conclusions were robust to the choice of interpolation or resampling method.

## Typical workflow

The code in this folder supports the following workflow:

```text
1. Initial coordinate alignment
   transcriptomic coordinates + metabolomic coordinates
   -> approximate transformation parameters
   -> initially transformed metabolomic coordinates

2. Fine image-based alignment
   initially aligned coordinates/images
   -> SURF-based feature matching and image-level refinement
   -> final registered metabolomic coordinates

3. Metabolite resampling
   final registered metabolomic coordinates + metabolite intensity matrix + transcriptomic bins
   -> metabolite values resampled onto transcriptomic spatial bins

4. Robustness checks
   area-overlap weighted resampling vs nearest-neighbour vs inverse-distance weighted resampling
   -> cluster-level and pairwise concordance analyses
```

## Example commands

From the repository root:

```bash
python code/run_initial_rigid_registration.py
python code/run_area_weighted_resampling.py
python code/run_nearest_idw_resampling.py
```

The example commands use small files in `example_data/` and write outputs to `outputs/`.

## Input and output notes

The initial alignment scripts require coordinate tables with at least `x` and `y` columns. The metabolite resampling scripts require registered MSI coordinates, transcriptomic spatial-bin coordinates and a metabolite intensity matrix. The exact input/output formats are summarized in `../docs/input_output_formats.md`.

The example files are provided only to demonstrate the expected data format and code execution. They are not the complete dataset used in the manuscript.

## References and external resources

- Bay, H., Tuytelaars, T. and Van Gool, L. SURF: Speeded Up Robust Features. European Conference on Computer Vision (ECCV), 2006. DOI: 10.1007/11744023_32. URL: http://www.vision.ee.ethz.ch/~surf/eccv06.pdf
- Bay, H., Ess, A., Tuytelaars, T. and Van Gool, L. Speeded-Up Robust Features (SURF). Computer Vision and Image Understanding 110, 346-359 (2008). DOI: 10.1016/j.cviu.2007.09.014. URL: https://doi.org/10.1016/j.cviu.2007.09.014
- SURF project page: http://www.vision.ee.ethz.ch/~surf
- OpenCV SURF tutorial: https://docs.opencv.org/4.x/df/dd2/tutorial_py_surf_intro.html
- SciPy Powell optimization documentation: https://docs.scipy.org/doc/scipy/reference/optimize.minimize-powell.html
- Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods 17, 261-272 (2020). DOI: 10.1038/s41592-019-0686-2.
