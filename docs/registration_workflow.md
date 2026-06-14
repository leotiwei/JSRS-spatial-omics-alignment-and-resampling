# Registration workflow

The joint spatial registration strategy (JSRS) starts from processed spatial transcriptomic coordinates and MALDI-MSI coordinates. The spatial transcriptomic coordinate system is used as the reference coordinate system.

## Workflow overview

1. **Coordinate preparation**
   - Extract RNA / Stereo-seq coordinates after transcriptomic processing.
   - Extract MALDI-MSI pixel coordinates and metabolite intensity matrices.
   - Standardize coordinate column names and section IDs.

2. **Initial rigid alignment**
   - Perform approximate rigid alignment using rotation and translation.
   - The example implementation in `code/jsrs/initial_rigid_registration.py` estimates an initial transformation by minimizing a rasterized point-cloud loss with `scipy.optimize.minimize` using the Powell algorithm.

3. **Image-based registration refinement**
   - Refine the initial alignment using SURF-based feature matching on section images or coordinate-derived reference images.
   - The refined transformation maps MALDI-MSI coordinates into the RNA coordinate system.

4. **Metabolite resampling**
   - Resample registered MSI intensities onto transcriptomic spatial units using the primary area-overlap weighted method.
   - Evaluate method robustness using nearest-neighbour and inverse-distance weighted resampling.

5. **Registration QC**
   - Quantify registration accuracy using predefined anatomical landmarks and tissue-mask overlap.
   - QC scripts are provided in `analysis/registration_qc`.

## Notes

The initial rigid alignment code is included to document the custom coordinate transformation step. The downstream refinement uses a standard SURF-based feature-matching workflow. Registration quality is reported using landmark-based error and mask-overlap metrics rather than relying only on visual overlays.
