# Aperture photometry Pipeline

Aperture photometry pipeline for ground-based and HST imaging data. The tool retrieves photometric reference catalogues from Vizier, builds a local comparison-star sequence, determines a photometric zeropoint, and measures calibrated aperture magnitudes.

Version 2026-07-28

---

## Contents

| Script | Purpose |
|--------|---------|
| `field_calibration.py` | Download reference catalogues from Vizier / NOIR DataLab |
| `improve_astrometry.py` | Improve WCS solution with astrometry.net |
| `align_images.py` | Align WCS of a science image to a reference image with SCAMP |
| `photometry.py` | Aperture photometry for ground-based images |
| `photometry_hst.py` | Aperture photometry + curve-of-growth for HST images |

| Module | Purpose |
|--------|---------|
| `extraction.py` | `sep` source detection, aperture photometry, forced photometry, background estimation |
| `calibration.py` | Local sequence, zeropoint bootstrap, science output tables, poststamp |
| `hst_routines.py` | HST-specific photometry, aperture corrections (PySynphot), curve of growth |
| `cat_tools.py` | Catalogue metadata, colour transformations, sky cross-matching, NOIR DataLab access |
| `fits_tools.py` | FITS utilities and WCS coordinate transformations |
| `utils.py` | Terminal colour codes, bootstrap statistics, and `setup_logging` |
| `plotsettings.py` | Matplotlib style |
| `sip_to_pv.py` | Convert SIP distortion keywords to PV format |

---

## Requirements

Install the `photometry` conda environment:

```bash
conda create -n photometry python=3.11
conda activate photometry
pip install numpy astropy scipy matplotlib sep photutils astroquery pysynphot
# For Legacy Survey queries (optional):
pip install noirlab-datalab
# For astrometry improvement (requires astrometry.net binary):
# macOS: brew install astrometry-net
conda install conda-forge::astrometry
# For image alignment (requires SCAMP binary):
conda install conda-forge::astromatic-scamp
```

**Key versions tested:** numpy 2.4.4, astropy 7.2.0, sep 1.4.1, photutils 3.0.0, pysynphot 2.0.0.

The pipeline no longer requires a SExtractor binary — source detection and photometry use the `sep` Python library (C-level SExtractor algorithms).

---

## Quick start

Run all commands from the directory containing the FITS files.

### Negative declinations

When the declination is negative, use the `=` form to prevent the shell interpreting the minus sign as a flag:

```bash
# Correct — works in bash, zsh, and sh
--dec="-12:20:20.2"
--dec="-0.726"

# Also works (leading space before the minus)
--dec " -12:20:20.2"
```

### 1. Retrieve reference catalogues

Downloads SDSS, PanSTARRS, DES, 2MASS, WISE catalogues and derives colour-transformed versions (Bessel, GROND, ZTF, HSC, DES systems).

```bash
python field_calibration.py --ra 11:33:41.550 --dec 00:43:33.50
```

Output files are written to `results/` (configurable with `--outdir`).
Catalogue naming: `{SOURCE}_{SYSTEM}_{FILTER}.cat`, e.g. `PS1_SDSS_r.cat`.

### 2. Improve astrometry (optional)

Refine the WCS solution using astrometry.net and convert SIP distortion terms to PV format for Source Extractor:

```bash
python improve_astrometry.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits IMAGE.fits
```

Output: `IMAGE_wcs.fits` with improved WCS solution.
Log: `IMAGE_astrometry.log` with WCS quality statistics (RMS residual, pixel scale, matched stars).

### 2b. Align images to a reference (optional)

Align the WCS of a science image to a reference image using SCAMP.  Useful when several epochs need a common astrometric reference frame:

```bash
python align_images.py \
    --ref-image reference.fits \
    --new-image science.fits
```

Output: `science_astro.fits` with the aligned WCS applied.
Log: `science_alignment.log` with SCAMP quality metrics.

The reference image must have a valid WCS.  Sources are detected with `sep` and windowed centroids are used for both the reference and input catalogs.

### 3. Get photometry
### 3a. Ground-based image

Requires a reference catalogue produced by `field_calibration.py`:

```bash
python photometry.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits SN2015bn_SDSS_r_wcs.fits \
    --ref-file results/SDSS_SDSS_r.cat
```

### 3b. HST data

```bash
python photometry_hst.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits SN2015bn_F625W_drc.fits
```

Requires pysynphot and a valid `$PYSYN_CDBS` environment variable pointing to the Calibration Data Base System (CDBS) directory.

---

## Photometry workflow

```
field_calibration.py
        | reference catalogue (.cat)
        v
photometry.py
  Step 1: Verify target in image footprint (WCS)
  Step 2: Build local comparison-star sequence
          +-- Run sep on reference star positions (ASSOC mode)
          +-- IQR-clip, flag filter, maxstars limit
          +-- Cross-match sep catalog and reference catalog
  Step 3: Zeropoint determination
          +-- Bootstrap resampling (NITER=1000, ~0.07 s)
          +-- Diagnostic plots: ZP = m_cat - m_inst vs. magnitude, FWHM distribution
  Step 4: Aperture photometry on science target
          +-- Run sep on full image (all sources)
          +-- Run sep with ASSOC near target position
          +-- If undetected: forced photometry at the expected position
              (same global background, RMS map, gain, and aperture radii
               as the detection path) -> 3-sigma upper limits
  Step 5: Calibrate and summarise
  Step 6: Poststamp image
  Step 7: Write output files
  Step 8: Clean up temporaries
```

---

## Output files

| File | Format | Description |
|------|--------|-------------|
| `*_astro.fits` | FITS | Aligned image with updated WCS (from `align_images.py`) |
| `*_alignment.log` | text | SCAMP quality metrics and run record (from `align_images.py`) |
| `*_phot.log` | ECSV | Calibrated science photometry summary |
| `*_zp.log` | ECSV | Zeropoint table (one row per aperture) |
| `*_all_abs_cal.fits` | FITS BINTABLE | Full source catalogue with calibrated magnitudes |
| `*_zp.pdf` | PDF | Zeropoint diagnostic plot (ZP = m_cat - m_inst vs. magnitude) |
| `*_fwhm.pdf` | PDF | FWHM distribution of comparison stars |
| `*_std.pdf` | PDF | Instrumental vs. apparent magnitude scatter plot |
| `*_poststamp.pdf` | PDF | Two-panel cutout: expected vs. observed position |

ECSV (Enhanced Character Separated Values) is plain text with an embedded YAML header that preserves column dtypes, units, and descriptions.  Files can be opened in any text editor and read back with `astropy.io.ascii.read()`. The FITS catalog is directly readable by TopCat, DS9, and `astropy.table.Table.read()`.

### Photometry log format

The `*_phot.log` ECSV file contains one row per measured quantity:

```
# %ECSV 1.0
# ---
# datatype: [{name: PROPERTY, datatype: string}, ...]
PROPERTY           VALUE     ERROR+   ERROR-   COMMENT
FILENAME           nan       nan      nan      IMAGE.fits
DATE-OBS           nan       nan      nan      2015-08-15T04:30:00
MJD                nan       nan      nan      57249.1875
EXPTIME            nan       nan      nan      300.0
RA                 173.423   nan      nan      degree
DEC                0.726     nan      nan      degree
XWIN_IMAGE_EXP     1212.0    nan      nan      pixel    <- WCS-expected position
YWIN_IMAGE_EXP     1211.0    nan      nan      pixel
XWIN_IMAGE_OBS     1212.0    nan      nan      pixel    <- detected centroid
YWIN_IMAGE_OBS     1211.0    nan      nan      pixel
DISTANCE (arcsec)  0.32      nan      nan      arcsec
MAG_AUTO           22.50     0.10     0.09     mag      <- Kron magnitude
MAG_PETRO          22.51     0.11     0.10     mag      <- Petrosian magnitude
MAG_APER           22.50     0.03     0.03     mag      <- smallest aperture
MAG_APER_1         22.30     0.02     0.02     mag
MAG_APER_2         22.23     0.02     0.02     mag
MAG_APER_3         21.75     0.02     0.02     mag      <- largest circular aperture
```

Whether the photometry is in the AB or Vega system depends on the source catalogue.

### Non-detections (forced photometry)

If no source is found within `--host-offset` of the target, the pipeline falls
back to **forced photometry** at the WCS-expected position instead of leaving
the row blank.  The forced measurement deliberately reuses the same recipe as
the detection path — global-mesh background subtraction, the per-pixel
background-RMS map as the error image, the same detector gain, and the same
aperture radii — so the reported brightness and uncertainty are directly
comparable with detected sources.  The `*_phot.log` then contains, per aperture:

```
MAG_APER_PHOTUTILS          26.66  nan  nan  mag      <- flux measurement (= 3-sigma
MAG_APER_PHOTUTILS_3SIGMA   26.66  nan  nan  mag         limit when flux <= 0)
FNU_APER_PHOTUTILS       -0.039  0.026  ...  microJy  <- flux density (may be negative)
```

For a non-detection the flux scatters around zero, so `MAG_APER_PHOTUTILS`
equals the 3-sigma upper limit `MAG_APER_PHOTUTILS_3SIGMA`, and the asymmetric
`ERROR+/ERROR-` columns are `nan`.  Smaller apertures give the deepest limits
(less sky noise integrated).

### Zeropoint table (`*_zp.log`)

One row per method (`MAG_AUTO`, `MAG_PETRO`, `MAG_APER[_N]`) with the zeropoint
`ZP` and its bootstrap errors, the aperture size (`diam(px)`, `diam(arcsec)`),
the aperture correction `AP_cor` (relative to the largest aperture), and
`MAG_3UL_GLOB` — a **global 3-sigma limiting magnitude** estimated as the median
calibrated magnitude of sources detected at ~3 sigma (`MAGERR` in 0.30-0.34).
It is an image-wide, catalogue-based depth indicator, complementary to the
sky-noise-based `*_3SIGMA` limit measured at the target position.

---

## Key options

### `field_calibration.py`

| Option | Default | Description |
|--------|---------|-------------|
| `--ra` | required | RA (HMS or decimal degrees) |
| `--dec` | required | Dec (DMS or decimal degrees) |
| `--radius` | 10 | Search radius (arcmin) |
| `--outdir` | `results/` | Output directory |

### `align_images.py`

| Option | Default | Description |
|--------|---------|-------------|
| `--ref-image` | required | Reference FITS image (must have a valid WCS) |
| `--new-image` | required | Science image to be aligned |
| `--det-thresh` | 5.0 | Source detection threshold (sigma above background) |
| `--distort-degrees` | 1 | SCAMP polynomial degree for the distortion model |
| `--position-maxerr` | 0.15 | Maximum positional offset for source matching (deg) |
| `--sn-thresholds` | `10.0,40.0` | SCAMP S/N thresholds for two-pass matching |
| `--keep-temp` | False | Keep intermediate LDAC catalogs and SCAMP files |

### `photometry.py`

| Option | Default | Description |
|--------|---------|-------------|
| `--ra` | required | RA of target |
| `--dec` | required | Dec of target |
| `--fits` | required | Science FITS file |
| `--ref-file` | required | Reference catalogue file (generate with field_calibration.py) |
| `--ap-diam` | `1.0 1.5 2.0 3.0` | Aperture diameters in units of FWHM |
| `--host-offset` | 5 | Detection radius around target (arcsec) |
| `--tol` | 1 | Cross-matching tolerance (arcsec) |
| `--gain` | auto | Extra header keyword to try for the CCD gain (e-/ADU). Even without this flag the standard keywords `CCDGAIN`, `ADCGAIN`, `ATODGAIN`, `GAIN`, `EGAIN` are searched automatically; if none is found the gain falls back to 1.0 with a yellow warning. |
| `--auto` | False | Non-interactive (auto magnitude cuts) |
| `--outdir` | `results/` | Output directory |

### `photometry_hst.py`

| Option | Default | Description |
|--------|---------|-------------|
| `--ra` | required | RA of target |
| `--dec` | required | Dec of target |
| `--fits` | required | HST drizzled FITS file |
| `--ap-diam` | `0.25 ... 3.00` | Aperture diameters in arcsec |
| `--ap-inner-annulus` | 1.25 | Inner background annulus (x aperture diameter) |
| `--ap-outer-annulus` | 2.5 | Outer background annulus (x aperture diameter) |
| `--centroid` | False | Recentre on nearest detected source |
| `--tol` | 2 | Centroid search radius (arcsec) |
| `--outdir` | `results/` | Output directory |

---

## Supported reference catalogues

| Catalogue | System | Filters | Star selection |
|-----------|--------|---------|----------------|
| SDSS DR12 | AB | u g r i z | class=6 (star), mode=1 (primary) |
| PanSTARRS DR1 | AB | g r i z y | Cross-matched with Gaia astrometric stars |
| DES DR1 | AB | g r i z Y | Morphological star classifier > 0.85 |
| SkyMapper DR2 | AB | u v g r i z | class_star > 0.95, flags=0 |
| 2MASS | Vega | J H K | --- |
| WISE (AllWISE) | Vega | W1 W2 | --- |
| Legacy Survey DR10 | AB | g r i z | PSF-type sources, S/N >= 3 |

### Colour transformations available

Starting from SDSS photometry (or equivalent via PS1->SDSS conversion):

- SDSS -> Bessel (B V R I)
- SDSS -> GROND (g r i z)
- SDSS -> DES (g r i z)
- PS1 -> SDSS (u g r i z)
- PS1 -> HSC (g r i z y)
- PS1 -> ZTF (g r)
- PS1 -> 2MASS J (synthetic; **stars only**)
- DES -> SDSS (g r i z)

**PS1 -> 2MASS J** is a synthetic transformation, $J = y_{P1} - P(g-i)$, derived by
synthetic photometry of the Pickles atlas through the PS1 and 2MASS passbands (see
[`ps1_to_2mass_J_transformation.ipynb`](ps1_to_2mass_J_transformation.ipynb)). It lets a
$J$-band zeropoint be set from deep PS1 photometry where 2MASS itself is too shallow. Only
valid for ordinary $\sim$B–K stars ($-0.86 < (g-i)_{P1} < 1.69$): the query drops M dwarfs
and later, very blue stars (colour range), and strong-metallicity / peculiar / binary /
bad-photometry outliers (a robust $(g-r, r-i)$ stellar-locus cut). Broadband PS1 colours
cannot separate luminosity class, so rare supergiants are not removed directly (they are
clipped in the zeropoint bootstrap). Written to `PS1_2MASS_J.cat`. Per-star scatter is
$\sim$0.05–0.1 mag, and the relation is robust to Milky-Way foreground reddening for
$E(B\!-\!V) < 0.3$.

---

## Performance

Performance on a 2423 x 2423 pixel SDSS r-band image:

| Step | Previous version (with sewpy/SExtractor) | New version (with sep) | Speedup |
|------|------------------------|-----------|---------|
| Source extraction (thresh=3-sigma) | ~30 s | **0.5 s** | ~60x |
| Catalog cross-match (10k x 10k) | ~2 s | **0.2 s** | ~10x |
| Bootstrap ZP (N=30 000) | ~30 s | **0.07 s** | ~430x |
| sky2xy (439 stars) | ~2 s | **0.02 s** | ~100x |

---

## Architecture notes

### Source detection (`sep`)

`sep` implements SExtractor's core algorithms as a Python/C library --- no external binary required.  Positions are 0-indexed internally and converted to 1-indexed FITS convention in the output table.

| SExtractor | sep equivalent |
|------------|---------------|
| XWIN_IMAGE | `objects['x'] + 1` |
| ALPHAWIN_J2000 | `wcs.wcs_pix2world(x, y, 0)` |
| MAG_AUTO | Kron aperture: `sep.kron_radius` + `sep.sum_ellipse` |
| MAG_PETRO | ~2x half-light radius: `sep.flux_radius` + `sep.sum_ellipse` |
| MAG_APER | Fixed circular aperture: `sep.sum_circle` |
| FLUX_RADIUS | `sep.flux_radius(frac=0.5)` --- 50% enclosed flux radius |
| FWHM_IMAGE | `2 x FLUX_RADIUS` (Gaussian approximation) |

### Cross-matching

All sky coordinate matching uses `astropy.coordinates.SkyCoord.match_to_catalog_sky()` which handles RA wrap-around correctly via a C-compiled KD-tree.

### Zeropoint determination

Bootstrap + Monte Carlo resampling with IQR outlier rejection (vectorised NumPy --- runs at ~430x the speed of the original Python loop). The diagnostic plot shows ZP = m_cat - m_inst (positive) vs. apparent
magnitude; outlier stars are highlighted in grey.

### Aperture photometry error model

Flux errors combine the sky-background noise and the source Poisson term in
quadrature, `fluxerr = sqrt(Sum_aperture rms^2 + flux/gain)`:

- **Background subtraction** uses `sep.Background` (SExtractor MESH): the image
  is tiled into `--back-size` (default 64 px) boxes, the background is
  sigma-clipped per box, the mesh is median-filtered (`--back-filtersize`), and
  spline-interpolated to full-resolution 2-D background and RMS maps.
- The **per-pixel RMS map** is passed to `sep.sum_circle`/`sum_ellipse` as the
  `err` image, so the sky term is included (not just source Poisson).
- The **gain** is read from the header (see `--gain`) so the Poisson term is
  correct; a missing gain falls back to 1.0 with a warning.
- Magnitude errors are asymmetric: `ERROR+/- = -/+2.5 log10((F -/+ dF)/F)`;
  when `dF > F` the bright-side error is `nan` (an effective non-detection).

**Caveat — correlated noise.** The aperture error assumes independent pixels.
Resampled/coadded images (e.g. SWarp output, non-integer effective gain) have
spatially correlated noise, so the true aperture error is somewhat larger than
`sqrt(Sum rms^2)`.  The HST path applies a drizzle correlation factor
(`corr_factor`); the ground-based `sep` path does not, matching SExtractor's
default behaviour.

### Forced photometry (non-detections)

`extraction.forced_photometry` measures at a fixed position with `sep` using the
identical background, RMS map, gain, and aperture radii as `extract_sources`, so
forced upper limits are consistent with detected-source magnitudes.  The
photutils-based `extraction.aperture_photometry` (local-annulus background plus a
drizzle `corr_factor`) is retained for the HST path (`photometry_hst.py`).

### Image alignment (`align_images.py`)

Source catalogs for SCAMP are built with `sep` (no SExtractor binary required).  Windowed centroids (`sep.winpos`) are used throughout.

The SCAMP quality log reports:

| Metric | Meaning |
|--------|---------|
| `XY_Contrast` | Pattern-match contrast; > 3 indicates a reliable, unambiguous solution |
| `dRA`, `dDec` | Mean residual positional offsets (arcsec) |
| `Sigma_alpha_int`, `Sigma_delta_int` | Internal 1-sigma scatter of matched pairs — primary precision indicator |
| `Chi2_int` | Reduced χ² of the internal fit; ≈ 1 is well-conditioned |

### Logging

Each script writes a log file (same stem as the input FITS, `.log` extension).  Log level is controlled via `--loglevel` (default: `INFO`). Python warnings (e.g. numpy RuntimeWarning) are routed into the same log file via `logging.captureWarnings`.  All configuration is handled by `utils.setup_logging`, which prevents duplicate log entries on module reload and sets `propagate=False`.

---

## Citation

If you use this pipeline, please cite the following software:

- **scamp**: [Arnouts (2006)](http://adsabs.harvard.edu/abs/2006ASPC..351..112B)
- **sep**: [Barbary (2016)](https://doi.org/10.21105/joss.00058)
- **astropy**: [Astropy Collaboration (2022)](https://doi.org/10.3847/1538-4357/ac7c74)
- **photutils**: [Bradley et al. (2022)](https://doi.org/10.5281/zenodo.6825092)


and 

```
@ARTICLE{Schulze2018a,
   author = {{Schulze}, S. and {Kr{\"u}hler}, T. and {Leloudas}, G. and {Gorosabel}, J. and 
	{Mehner}, A. and {Buchner}, J. and {Kim}, S. and {Ibar}, E. and 
	{Amor{\'{\i}}n}, R. and {Herrero-Illana}, R. and {Anderson}, J.~P. and 
	{Bauer}, F.~E. and {Christensen}, L. and {de Pasquale}, M. and 
	{de Ugarte Postigo}, A. and {Gallazzi}, A. and {Hjorth}, J. and 
	{Morrell}, N. and {Malesani}, D. and {Sparre}, M. and {Stalder}, B. and 
	{Stark}, A.~A. and {Th{\"o}ne}, C.~C. and {Wheeler}, J.~C.},
    title = "{Cosmic evolution and metal aversion in superluminous supernova host galaxies}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1612.05978},
 keywords = {galaxies: evolution, galaxies: high-redshift, galaxies: luminosity function, mass function, galaxies: starburst, galaxies: star formation},
     year = 2018,
    month = jan,
   volume = 473,
    pages = {1258},
      doi = {10.1093/mnras/stx2352},
   adsurl = {http://adsabs.harvard.edu/abs/2018MNRAS.473.1258S},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }
```

---

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Author

Steve Schulze (steve.schulze@weizmann.ac.il)