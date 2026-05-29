# Refactoring Plan — ENGRAVE Photometry Pipeline

**Date:** 2026-05-29  
**Baseline commit:** bf83153  
**Environment:** conda `photometry` (Python 3.11, numpy 2.4.4, astropy 7.2.0)

---

## Critical Issues (must fix for code to run)

### 1. sewpy → sep (BLOCKING)
`sewpy` is not installed. Replace with `sep 1.4.1` (already in env).

**Files affected:** `phot_routines.py`  
**Strategy:** Rewrite `sextractor_photometry()` using:
- `sep.Background()` for background estimation (replaces `photutils.make_source_mask` + `sigma_clipped_stats`)
- `sep.extract()` for source detection
- `sep.sum_circle()` for fixed-aperture photometry
- `sep.kron_radius()` + `sep.sum_ellipse()` for AUTO (Kron) photometry  
- `sep.flux_radius()` for half-light radius (→ FWHM)
- Post-extraction positional matching to replace ASSOC mode
- Save background-subtracted image as check image for poststamp

**Output table must preserve column names:**
`XWIN_IMAGE, YWIN_IMAGE, ALPHAWIN_J2000, DELTAWIN_J2000, FLAGS, A_IMAGE, B_IMAGE, THETA_IMAGE, FWHM_IMAGE, FWHM_WORLD, KRON_RADIUS, FLUX_RADIUS, FLUX_AUTO, FLUXERR_AUTO, FLUX_PETRO, FLUXERR_PETRO, MAG_AUTO, MAGERR_AUTO, MAG_PETRO, MAGERR_PETRO`  
Plus aperture columns: `MAG_APER, MAG_APER_1, ..., MAGERR_APER, ..., FLUX_APER, ..., FLUXERR_APER, ...`  
Plus ASSOC columns (when ASSOC_NAME provided): `VECTOR_ASSOC, NUMBER_ASSOC`

**Note on PHOT_APERTURES:** SExtractor convention — values are **diameters** in pixels. sep uses radii. Divide by 2.  
**Note on pixel indexing:** sep is 0-indexed; SExtractor/FITS is 1-indexed. Add 1 to x,y for output.  
**Note on ASSOC:** Input xy file is 1-indexed (from fits_tools.sky2xy). Subtract 1 when comparing to sep positions.

### 2. numpy 2.x `.view()` incompatibility (BLOCKING)
Pattern `np.asarray(table).view((float, n))` fails in numpy 2.x.

**Fix everywhere (3 locations):**
```python
# OLD (broken):
arr = np.asarray(table[keys]).view((float, len(table.dtype.names)))
# NEW:
from numpy.lib.recfunctions import structured_to_unstructured
arr = structured_to_unstructured(np.asarray(table[keys]), dtype=float)
```

**Locations:**
- `photometry.py` lines 499–500 and 614–615
- `phot_routines.py` `zeropoint()` lines 1523–1524

**Also fix `zeropoint()`:** Hardcodes exactly 4 aperture columns in `TABLE_NEW_keys`. Make dynamic:
```python
aper_keys = [k for k in TABLE_NEW.keys() if any(k.startswith(p) for p in 
             ['MAG_APER', 'MAGERR_APER', 'FLUX_APER', 'FLUXERR_APER'])]
base_keys = [k for k in required_base if k in TABLE_NEW.colnames]
TABLE_NEW_keys = base_keys + aper_keys
```

### 3. pysynphot — now installed (2.0.0)
`pysynphot` was missing but is now installed. No code changes needed for HST functions.

### 4. photutils 3.0 API changes
`photutils.make_source_mask` removed. Replace `background()` function with `sep.Background()`.  
Also fix: `photutils.CircularAperture` still works but import from `photutils.aperture` to be explicit.

---

## Bugs to Fix

### B1. EllipticalAperture overwrites CircularAperture
`phot_routines.py` lines 39–40: second line overwrites first. Remove EllipticalAperture line.

### B2. `sextractor_postprocess` indexing bug
`phot_routines.py` line 1493: After setting `DATA[key] <= 0` to NaN, the condition `DATA[key] <= 0` is always False (data is now NaN). Fix: save mask before modifying.

```python
# OLD:
DATA[key][DATA[key] <= 0.] = np.nan
DATA[key.replace('FLUX_', 'FLUXERR_')][DATA[key] <= 0] = np.nan  # BUG: key now has NaN
# NEW:
neg_mask = DATA[key] <= 0.
DATA[key][neg_mask] = np.nan
DATA[key.replace('FLUX_', 'FLUXERR_')][neg_mask] = np.nan
```

### B3. Print statement bug in `field_calibration.py`
Line 149: says "Convert SDSS -> GROND" but should say "Convert SDSS -> DES".

### B4. `import pdb` left in `fits_tools.py` line 1

### B5. Eccentricity bug in `photometry.py` (commented out)
Line 442: `np.sqrt(1 - (ref_stars['A_IMAGE'] / ref_stars['A_IMAGE'])**2)` — A/A is always 1. Should be B/A.

### B6. Gain keyword fallback missing
`phot_routines.py` `get_gain()` line 241–246: no error message before `sys.exit()` in the else branch. Add descriptive message.

### B7. `make_scicat` loop variable `i` incremented but not used
`phot_routines.py` line 1185: `i += 1` inside loop but loop variable is already managed by `enumerate()`. Remove.

---

## Performance Improvements

### P1. Cross-matching (biggest bottleneck for 10k-entry catalogs)
Replace custom KD-tree `wrapper_crossmatch()` with `astropy.coordinates.SkyCoord.match_to_catalog_sky()`.

```python
from astropy.coordinates import SkyCoord
import astropy.units as u

def crossmatch_sky(ra1, dec1, ra2, dec2, max_sep_arcsec):
    cat1 = SkyCoord(ra=ra1*u.deg, dec=dec1*u.deg)
    cat2 = SkyCoord(ra=ra2*u.deg, dec=dec2*u.deg)
    idx, d2d, _ = cat1.match_to_catalog_sky(cat2)
    mask = d2d.arcsec < max_sep_arcsec
    return idx, d2d.arcsec, mask
```

**Speedup:** astropy uses C-compiled KD-tree with proper RA wrap-around handling. For 10k entries, ~5-10x faster.

Keep `crossmatch_angular()` and `wrapper_crossmatch()` as thin wrappers for backward compatibility.

### P2. Vectorize `statNclip()` in `stat_tools.py`
Replace the Python while loop with numpy vectorized operations:
```python
# OLD: 1000-iteration Python loop (slow)
while i < NITER:
    random_MC = np.random.normal(A[mask,0], A[mask,1])
    ...

# NEW: vectorized
random_MC = np.random.normal(A[mask,0], A[mask,1], size=(NITER, len(mask)))
if len(A[:,0]) > 10:
    boot_idx = np.random.randint(0, len(mask), size=(NITER, len(mask)))
    random_array = np.median(random_MC[np.arange(NITER)[:,None], boot_idx], axis=1)
else:
    random_array = np.median(random_MC, axis=1)
```
**Speedup:** ~50-100x for NITER=1000, ~500-1000x for NITER=30000.

### P3. Vectorize `sky2xy()` in `fits_tools.py`
Replace the row-by-row loop with batch WCS transformation:
```python
# OLD: for loop over all catalog entries
for ii in range(len(cat)):
    _sky2xy = sky2xy_helper(FITS, ra[ii], dec[ii])

# NEW: batch WCS
ra_arr = np.array(cat[cat_keys[0]], dtype=float)
dec_arr = np.array(cat[cat_keys[1]], dtype=float)
x_arr, y_arr = w.wcs_world2pix(ra_arr, dec_arr, 1)
# Filter to image footprint
mask = (x_arr > 0) & (x_arr < naxis1) & (y_arr > 0) & (y_arr < naxis2)
return np.column_stack([x_arr[mask], y_arr[mask]])
```
**Speedup:** For 10k entries, O(n) → O(1) calls to WCS.

---

## Code Quality Improvements

### Q1. Remove Python 2 compatibility code
In `phot_routines.py` `local_sequence()`:
- Remove all `if sys.version_info[0] == 2: raw_input(...)` blocks
- Keep only `input(...)`

### Q2. Replace bare `except:` with specific exceptions
**Locations:** 
- `field_calibration.py` lines 166, 350, 456, 549, 622, 652, 675
- `phot_routines.py` lines 194, 254, 479, etc.
- Replace with `except Exception as e:` at minimum, log `e`

### Q3. Replace `os.system()` with `subprocess` or `pathlib`
- `field_calibration.py` line 578: `os.system('rm ...')` → `Path(...).unlink()`
- `photometry.py` lines 924–931: use `pathlib.Path.unlink()` or `glob`
- `improve_astrometry.py` lines 356–359: use `subprocess.run()` or `pathlib`
- `cat_tools.py` line 332: `os.system('mkdir %s')` → `Path(OUTDIR).mkdir()`

### Q4. Remove dead code
- `phot_routines.py` lines 39: commented-out `import pdb`
- Numerous commented-out code blocks throughout
- Python 2 compatibility blocks

### Q5. Add type hints to all public functions

### Q6. Consistent logging — use `LOGGER.info()` instead of `print()` for progress messages

---

## Module Restructuring

### M1. Move `routines_noir.py` content into `cat_tools.py`
`routines_noir.py` is just another catalog query; consolidate into `cat_tools.py`.
Handle missing `dl` module gracefully with a try/except at import time.

### M2. `misc.py` — keep as is (tiny, used widely)

### M3. Add `__all__` to `plotsettings.py` to prevent star-import namespace pollution

---

## Diagnostic Plot Improvements

### D1. Add residual scatter plot to zeropoint diagnostics
Show (cat_mag - zp_calibrated_mag) vs. cat_mag to diagnose systematic trends.

### D2. Improve poststamp — add compass (N/E orientation)

### D3. Add a combined summary plot showing all aperture results

---

## File-by-File Summary of Changes

| File | Changes |
|------|---------|
| `phot_routines.py` | Replace sewpy→sep, fix numpy 2.x in `zeropoint()`, fix aperture bug, fix postprocess bug, remove Py2 compat, fix `background()` |
| `photometry.py` | Fix `.view()` numpy 2.x (lines 499–500, 614–615), `os.system`→`pathlib`, fix eccentricity comment |
| `fits_tools.py` | Remove `import pdb`, vectorize `sky2xy()` |
| `cat_tools.py` | Add astropy-based `crossmatch_sky()`, docstrings, `os.system`→`pathlib` |
| `stat_tools.py` | Vectorize `statNclip()` |
| `field_calibration.py` | Fix print bug (DES msg), fix bare excepts, `os.system`→`pathlib` |
| `improve_astrometry.py` | `os.system`→`pathlib`, improve error handling |
| `plotsettings.py` | Add `__all__` |
| `routines_noir.py` | Graceful `dl` import |
| `misc.py` | No changes |
| `README.md` | Full update with dependencies, examples, output description |

---

## Testing

Run from `test/` subdirectory:

```bash
# Field calibration (downloads catalog)
python ../field_calibration.py --ra 11:33:41.550 --dec 00:43:33.50

# Astrometry improvement
python ../improve_astrometry.py --ra 11:33:41.550 --dec 00:43:33.50 --fits SN2015bn_SDSS_r_wcs.fits

# Photometry (requires results/SDSS_SDSS_r.cat from field_calibration)
python ../photometry.py --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits SN2015bn_SDSS_r_wcs.fits \
    --ref-file results/SDSS_SDSS_r.cat

# HST photometry
python ../photometry_hst.py --ra 11:33:41.550 --dec 00:43:33.50 \
    --fits SN2015bn_F625W_drc.fits
```

---

## Profiling

Add cProfile wrapper to identify bottlenecks:
```python
# At top of photometry.py (or separate profile_photometry.py):
import cProfile, pstats
pr = cProfile.Profile()
pr.enable()
# ... main code ...
pr.disable()
stats = pstats.Stats(pr).sort_stats('cumulative')
stats.print_stats(20)
```

Expected bottlenecks (in order):
1. SExtractor/sep source extraction
2. Bootstrap ZP calculation (`statNclip` with NITER=30000)
3. Catalog cross-matching (10k entries)
4. VizieR queries (I/O bound)
