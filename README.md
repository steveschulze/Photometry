# ENGRAVE Photometry Pipeline

Aperture photometry pipeline for ground-based and HST imaging of transients
(supernovae, gamma-ray bursts, etc.).  Retrieves photometric reference
catalogues from Vizier, builds a local comparison-star sequence, determines
a robust photometric zeropoint, and measures calibrated aperture magnitudes.

Version 2026-05-29

---

## Contents

| Script | Purpose |
|--------|---------|
| `field_calibration.py` | Download reference catalogues from Vizier / NOIR DataLab |
| `improve_astrometry.py` | Improve WCS solution with astrometry.net |
| `photometry.py` | Aperture photometry for ground-based images |
| `photometry_hst.py` | Aperture photometry + curve-of-growth for HST images |


| Module | Purpose |
|--------|---------|
| `extraction.py` | sep source detection, aperture photometry, background estimation |
| `calibration.py` | Local sequence, zeropoint bootstrap, science output tables, poststamp |
| `hst_routines.py` | HST-specific photometry, aperture corrections (pysynphot), curve-of-growth |
| `cat_tools.py` | Catalogue metadata, colour transformations, sky cross-matching, NOIR DataLab |
| `fits_tools.py` | FITS utilities and WCS coordinate transformations |
| `utils.py` | Terminal colour codes and bootstrap statistics |
| `plotsettings.py` | Matplotlib style and Vigit colour scheme |
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
```

**Key versions tested:** numpy 2.4.4, astropy 7.2.0, sep 1.4.1, photutils 3.0.0, pysynphot 2.0.0.

The pipeline no longer requires a SExtractor binary — source detection and
photometry use the `sep` Python library (C-level SExtractor algorithms).

---

## Quick start

Run all commands from the directory containing the FITS files.

### 1. Retrieve reference catalogues

Downloads SDSS, PanSTARRS, DES, 2MASS, WISE catalogues and derives colour-
transformed versions (Bessel, GROND, ZTF, HSC, DES systems).

```bash
python field_calibration.py --ra 11:33:41.550 --dec 00:43:33.50
```

Output files are written to `results/` (configurable with `--outdir`).
Catalogue naming: `{SOURCE}_{SYSTEM}_{FILTER}.cat`, e.g. `PS1_SDSS_r.cat`.

### 2. Improve astrometry (optional)

Refine the WCS solution using astrometry.net and convert SIP distortion
terms to PV format:

```bash
python improve_astrometry.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits IMAGE.fits
```

Output: `IMAGE_wcs.fits` with improved WCS solution.

### 3. Run photometry

Requires a reference catalogue produced by `field_calibration.py`:

```bash
python photometry.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits SN2015bn_SDSS_r_wcs.fits \
    --ref-file results/SDSS_SDSS_r.cat
```

Or specify the catalogue by name (downloaded on-the-fly from Vizier):

```bash
python photometry.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits IMAGE.fits \
    --ref-cat SDSS --ref-filter r
```

### 4. HST photometry

```bash
python photometry_hst.py \
    --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits SN2015bn_F625W_drc.fits
```

Requires pysynphot and a valid `$PYSYN_CDBS` environment variable pointing
to the Calibration Data Base System (CDBS) directory.

---

## Photometry pipeline workflow

```
field_calibration.py
        │ reference catalogue (.cat)
        ▼
photometry.py
  Step 1: Verify target in image footprint (WCS)
  Step 2: Build local comparison-star sequence
          ├─ Run sep on reference star positions (ASSOC mode)
          ├─ IQR-clip, flag filter, maxstars limit
          └─ Cross-match sep catalog × reference catalog
  Step 3: Zeropoint determination
          ├─ Bootstrap resampling (NITER=1000, ~0.07 s)
          └─ Diagnostic plots: ZP vs. magnitude, FWHM distribution
  Step 4: Aperture photometry on science target
          ├─ Run sep on full image (all sources)
          └─ Run sep with ASSOC near target position
  Step 5: Calibrate and summarise
  Step 6: Poststamp image
  Step 7: Write output files
  Step 8: Clean up temporaries
```

---

## Output files

| File | Description |
|------|-------------|
| `*_phot.log` | Calibrated science photometry summary |
| `*_zp.log` | Zeropoint table (one row per aperture) |
| `*_all_abs_cal.phot` | Full source catalogue with calibrated magnitudes |
| `*_zp.pdf` | Zeropoint diagnostic plot (ZP vs. reference magnitude) |
| `*_fwhm.pdf` | FWHM distribution of comparison stars |
| `*_std.pdf` | Instrumental vs. apparent magnitude scatter plot |
| `*_poststamp.pdf` | Two-panel cutout: expected vs. observed position |

### Photometry log format

The `*_phot.log` file contains one row per measured quantity:

```
PROPERTY          VALUE     ERROR+   ERROR-   COMMENT
FILENAME          nan       nan      nan      IMAGE.fits
DATE-OBS          nan       nan      nan      2015-08-15T04:30:00
MJD               nan       nan      nan      57249.1875
EXPTIME           nan       nan      nan      300.0
RA                173.423   nan      nan      degree
DEC               0.726     nan      nan      degree
XWIN_IMAGE_OBS    1212.0    nan      nan      pixel
YWIN_IMAGE_OBS    1211.0    nan      nan      pixel
DISTANCE (arcsec) 0.32      nan      nan      arcsec
MAG_AUTO          22.50     0.10     0.09     mag
MAG_APER          22.50     0.03     0.03     mag    ← smallest aperture
MAG_APER_1        22.30     0.02     0.02     mag
MAG_APER_2        22.23     0.02     0.02     mag
MAG_APER_3        21.75     0.02     0.02     mag    ← largest aperture
```

All magnitudes are in the AB system.

---

## Key options

### `field_calibration.py`

| Option | Default | Description |
|--------|---------|-------------|
| `--ra` | required | RA (HMS or decimal degrees) |
| `--dec` | required | Dec (DMS or decimal degrees) |
| `--radius` | 10 | Search radius (arcmin) |
| `--outdir` | `results/` | Output directory |

### `photometry.py`

| Option | Default | Description |
|--------|---------|-------------|
| `--ra` | required | RA of target |
| `--dec` | required | Dec of target |
| `--fits` | required | Science FITS file |
| `--ref-file` | — | Reference catalogue file |
| `--ref-cat` | — | Catalogue name (SDSS, 2MASS, PS1) |
| `--ref-filter` | — | Filter of reference catalogue |
| `--ap-diam` | `1.0 1.5 2.0 3.0` | Aperture diameters in units of FWHM |
| `--host-offset` | 5 | Detection radius around target (arcsec) |
| `--tol` | 1 | Cross-matching tolerance (arcsec) |
| `--gain` | — | Header keyword for CCD gain (e⁻/ADU) |
| `--auto` | False | Non-interactive (auto magnitude cuts) |
| `--outdir` | `results/` | Output directory |

---

## Supported reference catalogues

| Catalogue | System | Filters | Star selection |
|-----------|--------|---------|----------------|
| SDSS DR12 | AB | u g r i z | class=6 (star), mode=1 (primary) |
| PanSTARRS DR1 | AB | g r i z y | Cross-matched with Gaia astrometric stars |
| DES DR1 | AB | g r i z Y | Morphological star classifier > 0.85 |
| SkyMapper DR2 | AB | u v g r i z | class_star > 0.95, flags=0 |
| 2MASS | Vega | J H K | — |
| WISE (AllWISE) | Vega | W1 W2 | — |
| Legacy Survey DR10 | AB | g r i z | PSF-type sources, S/N ≥ 3 |

### Colour transformations available

Starting from SDSS photometry (or equivalent via PS1→SDSS conversion):

- SDSS → Bessel (B V R I)
- SDSS → GROND (g r i z)
- SDSS → DES (g r i z)
- PS1 → SDSS (u g r i z)
- PS1 → HSC (g r i z y)
- PS1 → ZTF (g r)
- DES → SDSS (g r i z)

---

## Performance

Performance on a 2423 × 2423 pixel SDSS r-band image:

| Step | Old (sewpy/SExtractor) | New (sep) | Speedup |
|------|------------------------|-----------|---------|
| Source extraction (thresh=3σ) | ~30 s | **0.5 s** | ~60× |
| Catalog cross-match (10k × 10k) | ~2 s | **0.2 s** | ~10× |
| Bootstrap ZP (N=30 000) | ~30 s | **0.07 s** | ~430× |
| sky2xy (439 stars) | ~2 s | **0.02 s** | ~100× |

---

## Architecture notes

### Source detection (`sep`)

`sep` implements SExtractor's core algorithms as a Python/C library —
no external binary required.  Positions are 0-indexed internally and
converted to 1-indexed FITS convention in the output table.

| SExtractor | sep equivalent |
|------------|---------------|
| XWIN_IMAGE | `objects['x'] + 1` |
| ALPHAWIN_J2000 | `wcs.wcs_pix2world(x, y, 0)` |
| MAG_AUTO | Kron aperture: `sep.kron_radius` + `sep.sum_ellipse` |
| MAG_PETRO | ~2× half-light radius: `sep.flux_radius` + `sep.sum_ellipse` |
| MAG_APER | Fixed circular aperture: `sep.sum_circle` |
| FLUX_RADIUS | `sep.flux_radius(frac=0.5)` — 50% enclosed flux radius |
| FWHM_IMAGE | `2 × FLUX_RADIUS` (Gaussian approximation) |

### Cross-matching

All sky coordinate matching uses `astropy.coordinates.SkyCoord.match_to_catalog_sky()`
which handles RA wrap-around correctly via a C-compiled KD-tree.

### Zeropoint determination

Bootstrap + Monte Carlo resampling with IQR outlier rejection
(vectorised NumPy — runs at ~430× the speed of the original Python loop).

---

## Bugs fixed

- **EllipticalAperture overwrite**: forced photometry was using elliptical
  apertures instead of circular ones.
- **`sextractor_postprocess` mask bug**: negative-flux mask was applied after
  setting values to NaN, causing the error-column update to be silently skipped.
- **`field_calibration.py` DES print**: displayed "Convert SDSS → GROND" for
  the DES conversion step.
- **numpy 2.x `.view()` incompatibility**: replaced with
  `numpy.lib.recfunctions.structured_to_unstructured()`.
- **`photutils.make_source_mask` removed in photutils 3.0**: replaced with
  `sep.Background`.
- **Python 2 `raw_input()` remnants**: removed.
- **Bare `except:` clauses**: replaced with `except Exception as exc:`.
- **`os.system()` shell commands**: replaced with `pathlib` / `subprocess`.

---

## Citation

If you use this pipeline, please cite the following software:

- **sep**: [Barbary (2016)](https://doi.org/10.21105/joss.00058)
- **astropy**: [Astropy Collaboration (2022)](https://doi.org/10.3847/1538-4357/ac7c74)
- **photutils**: [Bradley et al. (2022)](https://doi.org/10.5281/zenodo.6825092)

---

## Author

Steve Schulze (steve.schulze@weizmann.ac.il)
