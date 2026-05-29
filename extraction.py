"""
extraction.py — Source detection and raw aperture photometry.

Wraps the sep library (SExtractor algorithms in Python/C) for background
estimation, source detection, and aperture/Kron/Petrosian photometry.
No external binary is required.

Public API
----------
background(image, ...)           Global background RMS via sep.Background
background_local(image, aper)    Sigma-clipped stats inside an annulus
get_gain(fits_path, keyword, ...) Read detector gain from FITS header
aperture_photometry(...)          Circular aperture photometry via photutils
extract_sources(fits_path, ...)   sep-based source detection + photometry
postprocess_catalog(catalog, ...) Convert negative fluxes; recompute errors
setup_sextractor(output_dir)      Write SExtractor config files for astrometry.net
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"


import logging
import sys
import warnings
from pathlib import Path
from typing import Sequence

import numpy as np
import sep
from astropy import stats, table, wcs
from astropy.io import ascii, fits
from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from photutils.aperture import (
    CircularAnnulus,
    CircularAperture,
    aperture_photometry as phot_aperture_photometry,
)

from utils import bcolors

warnings.filterwarnings("ignore", category=AstropyUserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Module-level SExtractor defaults (match SExtractor 2.x defaults)
# ---------------------------------------------------------------------------

KRON_FACTOR: float = 2.5    # PHOT_AUTOPARAMS[0]
KRON_MIN_R:  float = 3.5    # PHOT_AUTOPARAMS[1]  (minimum Kron aperture radius)
PETRO_FACTOR: float = 2.0   # PHOT_PETROPARAMS[0]
PETRO_MIN_R:  float = 3.5   # PHOT_PETROPARAMS[1]


# ---------------------------------------------------------------------------
# Background estimation
# ---------------------------------------------------------------------------

def background(
    image: np.ndarray,
    back_size: int = 64,
    back_filtersize: int = 3,
) -> float:
    """Estimate the global background RMS of an image using sep.

    Parameters
    ----------
    image : np.ndarray
        2-D science image array.
    back_size : int
        Background mesh size in pixels (default: 64).
    back_filtersize : int
        Background filter size (default: 3).

    Returns
    -------
    float
        Global background RMS.
    """
    data = np.ascontiguousarray(image.astype(np.float64))
    bkg  = sep.Background(data, bw=back_size, bh=back_size,
                          fw=back_filtersize, fh=back_filtersize)
    return float(bkg.globalrms)


def background_local(
    image: np.ndarray,
    aperture: CircularAnnulus,
) -> tuple[float, float, float]:
    """Compute sigma-clipped statistics inside a background annulus.

    Parameters
    ----------
    image : np.ndarray
        2-D science image array.
    aperture : photutils.aperture.CircularAnnulus
        Background annulus aperture object.

    Returns
    -------
    tuple[float, float, float]
        (mean, median, std) of pixels inside the aperture.
    """
    bkg_mask = aperture.to_mask(method='center')
    bkg_data = bkg_mask.multiply(image)
    aper_data = bkg_data[bkg_mask.data > 0]
    return tuple(stats.sigma_clipped_stats(aper_data))  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Gain extraction
# ---------------------------------------------------------------------------

def get_gain(
    fits_path: str | Path,
    keyword: str | None,
    logger: logging.Logger | None = None,
) -> float:
    """Read detector gain from the FITS header.

    Searches for common gain keywords in priority order:
    CCDGAIN, ADCGAIN, ATODGAIN, then the user-supplied *keyword*.

    Parameters
    ----------
    fits_path : str | Path
        Path to the FITS file.
    keyword : str | None
        Header keyword to try last. Pass ``None`` to skip and return 1.
    logger : logging.Logger | None
        Logger for error reporting.

    Returns
    -------
    float
        Gain in e⁻/ADU, or 1.0 if *keyword* is None.
    """
    if keyword is None:
        return 1.0

    try:
        with fits.open(str(fits_path)) as hdu:
            header = hdu[0].header.copy()
            if len(hdu) > 1:
                header += hdu[1].header
            if len(hdu) > 2:
                header += hdu[2].header

        for key in ('CCDGAIN', 'ADCGAIN', 'ATODGAIN', keyword):
            if key in header:
                return float(header[key])

        msg = f'Gain keyword ({keyword}) not found in {fits_path}'
        print(bcolors.BOLD + bcolors.FAIL + f'\n{msg}\n' + bcolors.ENDC)
        if logger:
            logger.error(msg)
        sys.exit(1)

    except (OSError, KeyError) as exc:
        msg = f'Failed to read gain from {fits_path}: {exc}'
        print(bcolors.BOLD + bcolors.FAIL + f'\n{msg}\n' + bcolors.ENDC)
        if logger:
            logger.error(msg)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Aperture photometry (photutils-based, used for forced photometry)
# ---------------------------------------------------------------------------

def aperture_photometry(
    image: np.ndarray,
    positions: Sequence[float] | np.ndarray,
    radii: np.ndarray,
    inner_annulus: np.ndarray,
    outer_annulus: np.ndarray,
    rms: float,
    gain: float = 1.0,
    corr_factor: float = 1.0,
    zeropoint: float | np.ndarray = 0.0,
) -> table.Table:
    """Circular aperture photometry for one target at multiple aperture sizes.

    Subtracts a local background measured in an annulus around each aperture,
    computes shot and background noise separately, and returns magnitudes in
    the AB system plus flux densities in micro-Jy.

    Pixel positions follow the 1-indexed FITS convention (photutils uses
    0-indexed internally; the conversion is applied here).

    Parameters
    ----------
    image : np.ndarray
        Background-subtracted 2-D image array.
    positions : array-like, shape (2,)
        (x, y) pixel position of the target (1-indexed).
    radii : np.ndarray
        Aperture radii in pixels.
    inner_annulus : np.ndarray
        Inner background annulus radii (pixels).
    outer_annulus : np.ndarray
        Outer background annulus radii (pixels).
    rms : float
        Global image RMS (used as a fallback; local RMS is preferred).
    gain : float
        Effective detector gain in e⁻/ADU (default: 1).
    corr_factor : float
        Correlated-noise correction factor FA (default: 1 = uncorrelated).
    zeropoint : float or np.ndarray
        AB magnitude zeropoint per aperture.

    Returns
    -------
    astropy.table.Table
        One row per target with columns for each aperture:
        ``bkg_sub_aperture_sum_N``, ``TOTAL_ERROR_N``, ``FNU_APER_N``,
        ``FNUERR_APER_N``, ``MAG_APER_N``, ``MAGERRP_APER_N``,
        ``MAGERRM_APER_N``, ``ZP_APER_N``.
    """
    # photutils uses 0-indexed positions
    pos_0 = np.asarray(positions, dtype=float) - 1.0
    zp    = np.broadcast_to(np.asarray(zeropoint, dtype=float), len(radii))

    src_apers = [CircularAperture(pos_0, r=r) for r in radii]
    bkg_apers = [CircularAnnulus(pos_0,
                                  r_in=inner_annulus[i],
                                  r_out=outer_annulus[i])
                 for i in range(len(radii))]

    src_tbl = phot_aperture_photometry(image, src_apers)
    bkg_tbl = phot_aperture_photometry(image, bkg_apers)

    if 'aperture_sum_0' not in src_tbl.colnames:
        src_tbl.rename_column('aperture_sum', 'aperture_sum_0')
        bkg_tbl.rename_column('aperture_sum', 'aperture_sum_0')

    local_rms   = [background_local(image, bkg_apers[i])[-1]
                   for i in range(len(radii))]
    area_ratios = [src_apers[i].area / bkg_apers[i].area
                   for i in range(len(radii))]

    for i, (ratio, lrms, zp_i) in enumerate(zip(area_ratios, local_rms, zp)):
        bkg_scaled = bkg_tbl[f'aperture_sum_{i}'] * ratio
        net_flux   = src_tbl[f'aperture_sum_{i}'] - bkg_scaled

        src_tbl[f'bkg_{i}']                    = bkg_scaled
        src_tbl[f'bkg_sub_aperture_sum_{i}']   = net_flux
        src_tbl[f'ZP_APER_{i}']                = float(zp_i)

        shot_noise = np.sqrt(np.abs(net_flux) / gain)
        bkg_noise  = np.sqrt(lrms**2 * src_apers[i].area / corr_factor)

        src_tbl[f'SHOT_NOISE_{i}']  = shot_noise
        src_tbl[f'BKG_NOISE_{i}']   = bkg_noise
        src_tbl[f'TOTAL_ERROR_{i}'] = np.sqrt(shot_noise**2 + bkg_noise**2)

        # Flux density in micro-Jy (AB ZP 23.9 ↔ 1 µJy)
        factor   = 10.0**(-0.4 * (float(zp_i) - 23.9))
        fnu      = net_flux * factor
        fnu_err  = src_tbl[f'TOTAL_ERROR_{i}'] * factor

        src_tbl[f'FNU_APER_{i}']    = fnu
        src_tbl[f'FNUERR_APER_{i}'] = fnu_err

        with np.errstate(divide='ignore', invalid='ignore'):
            mag = np.where(
                fnu > 0,
                -2.5 * np.log10(np.maximum(fnu, 1e-30)) + 23.9,
                -2.5 * np.log10(np.maximum(3 * fnu_err, 1e-30)) + 23.9)
            errp = np.where(fnu > 0,
                            -2.5 * np.log10(fnu - fnu_err)
                            + 2.5 * np.log10(fnu), np.nan)
            errm = np.where(fnu > 0,
                            -2.5 * np.log10(fnu)
                            + 2.5 * np.log10(fnu + fnu_err), np.nan)

        src_tbl[f'MAG_APER_{i}']     = mag
        src_tbl[f'MAGERRP_APER_{i}'] = errp
        src_tbl[f'MAGERRM_APER_{i}'] = errm

    for col in src_tbl.colnames:
        if 'MAG_APER' in col or 'aperture' in col:
            src_tbl[col].info.format = '%.3f'
        elif 'FNU' in col or 'NOISE' in col or 'TOTAL' in col:
            src_tbl[col].info.format = '%.3e'

    return src_tbl


# ---------------------------------------------------------------------------
# FITS loader helper
# ---------------------------------------------------------------------------

def _load_fits_data(fits_path: str | Path) -> tuple[np.ndarray, fits.Header]:
    """Load image data and combined header from a FITS file.

    Handles single-extension and multi-extension FITS (e.g. HST DRZ).
    Returns data as a native-byte-order, C-contiguous float64 array.

    Parameters
    ----------
    fits_path : str | Path
        Path to the FITS file.

    Returns
    -------
    tuple[np.ndarray, fits.Header]
        ``(data, header)`` — image array and combined primary+extension header.
    """
    with fits.open(str(fits_path)) as hdulist:
        header = hdulist[0].header.copy()
        if len(hdulist) > 1:
            try:
                _ = hdulist[1].data.shape  # test data presence
                header += hdulist[1].header
                data    = hdulist[1].data
            except (AttributeError, TypeError):
                data = hdulist[0].data
        else:
            data = hdulist[0].data

    data = np.array(data, dtype=np.float64)
    if not data.flags['C_CONTIGUOUS']:
        data = np.ascontiguousarray(data)
    return data, header


# ---------------------------------------------------------------------------
# SExtractor config writer (for astrometry.net compatibility)
# ---------------------------------------------------------------------------

def setup_sextractor(output_dir: str | Path = '.') -> None:
    """Write default SExtractor configuration files.

    Writes ``default.conv``, ``default.nnw``, ``default.param``, and
    ``default.sex`` into *output_dir*.  These files are consumed by
    astrometry.net's ``solve-field`` command.

    Parameters
    ----------
    output_dir : str | Path
        Target directory (default: current directory).
    """
    d = Path(output_dir)

    files: dict[str, str] = {
        'default.conv': (
            "CONV NORM\n"
            "# 3x3 all-ground convolution mask, FWHM = 2 pixels.\n"
            "1 2 1\n2 4 2\n1 2 1\n"
        ),
        'default.nnw': (
            "NNW\n"
            " 3 10 10  1\n"
            "-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01"
            " -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01\n"
            " 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00"
            "  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00\n"
            "-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00"
            "  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01"
            " -2.32028e+00\n"
            " 0.00000e+00\n 1.00000e+00\n"
        ),
        'default.param': (
            "NUMBER\nXWIN_IMAGE\nYWIN_IMAGE\n"
            "ALPHAWIN_J2000\nDELTAWIN_J2000\n"
            "FLAGS\nFWHM_IMAGE\nCLASS_STAR\n"
            "FLUX_AUTO\nFLUXERR_AUTO\nFLUX_RADIUS\n"
            "A_IMAGE\nB_IMAGE\nTHETA_IMAGE\nKRON_RADIUS\n"
        ),
        'default.sex': (
            "CATALOG_NAME     test.cat\n"
            "CATALOG_TYPE     ASCII_HEAD\n"
            "PARAMETERS_NAME  default.param\n"
            "DETECT_TYPE      CCD\n"
            "DETECT_MINAREA   3\n"
            "DETECT_THRESH    1.5\n"
            "ANALYSIS_THRESH  1.5\n"
            "FILTER           Y\n"
            "FILTER_NAME      default.conv\n"
            "DEBLEND_NTHRESH  32\n"
            "DEBLEND_MINCONT  0.005\n"
            "CLEAN            Y\n"
            "CLEAN_PARAM      1.0\n"
            "WEIGHT_TYPE      NONE\n"
            "PHOT_APERTURES   5\n"
            "PHOT_AUTOPARAMS  2.5, 3.5\n"
            "PHOT_PETROPARAMS 2.0, 3.5\n"
            "PHOT_FLUXFRAC    0.5\n"
            "SATUR_LEVEL      50000.0\n"
            "SATUR_KEY        SATURATE\n"
            "MAG_ZEROPOINT    0.0\n"
            "GAIN             0.0\n"
            "GAIN_KEY         GAIN\n"
            "PIXEL_SCALE      1.0\n"
            "SEEING_FWHM      1.2\n"
            "STARNNW_NAME     default.nnw\n"
            "BACK_TYPE        AUTO\n"
            "BACK_SIZE        64\n"
            "BACK_FILTERSIZE  3\n"
            "CHECKIMAGE_TYPE  NONE\n"
            "MEMORY_OBJSTACK  3000\n"
            "MEMORY_PIXSTACK  300000\n"
            "MEMORY_BUFSIZE   1024\n"
            "VERBOSE_TYPE     QUIET\n"
            "WRITE_XML        N\n"
        ),
    }

    for fname, content in files.items():
        (d / fname).write_text(content)


# ---------------------------------------------------------------------------
# Main source extraction + photometry
# ---------------------------------------------------------------------------

def extract_sources(
    fits_path: str | Path,
    *,
    analysis_thresh: float = 1.0,
    assoc_file: str | Path | None = None,
    assoc_cols: str = "1,2",
    assoc_radius: float = 10.0,
    back_size: int = 64,
    back_filtersize: int = 3,
    deblend_nthresh: int = 32,
    deblend_mincont: float = 0.005,
    detect_thresh: float = 1.0,
    flag: str = "",
    gain_key: str | None = None,
    logger: logging.Logger | None = None,
    output_dir: str | Path | None = None,
    phot_apertures: np.ndarray | None = None,
    ref_image: str | Path = "",
) -> table.Table:
    """Detect sources and measure photometry using sep.

    Replaces the original sewpy/SExtractor backend.  The returned table
    uses the same column names as the original SExtractor FITS_LDAC output
    so that all downstream code is unaffected.

    Coordinates follow the 1-indexed FITS/SExtractor convention in all
    output columns (XWIN_IMAGE, YWIN_IMAGE).  Internally sep uses 0-indexed
    positions; the conversion is applied before return.

    Detection threshold
    -------------------
    Uses a global-RMS-based absolute threshold (thresh = N_sigma × globalrms)
    rather than passing a per-pixel error map to sep.  This matches
    SExtractor's default behaviour (BACK_TYPE=AUTO, no weight map).

    ASSOC mode
    ----------
    When *assoc_file* is supplied, only sources within *assoc_radius* pixels
    of any position in that file are returned — the same semantics as
    SExtractor ASSOCSELEC_TYPE=MATCHED.

    Aperture columns
    ----------------
    The first aperture has no numeric suffix (``MAG_APER``, ``FLUX_APER``);
    subsequent apertures are indexed from 1 (``MAG_APER_1``, …).

    Parameters
    ----------
    fits_path : str | Path
        Science FITS file.
    analysis_thresh : float
        Analysis threshold in sigma (default: 1.0).
    assoc_file : str | Path | None
        ASCII file with 1-indexed pixel (x, y) positions for cross-matching.
    assoc_cols : str
        Column numbers (1-based) for x, y [, value] in *assoc_file* (default: "1,2").
    assoc_radius : float
        Matching radius in pixels (default: 10.0).
    back_size : int
        Background mesh size (default: 64 pixels).
    back_filtersize : int
        Background filter size (default: 3).
    deblend_nthresh : int
        Deblending sub-thresholds (default: 32; SExtractor default).
    deblend_mincont : float
        Minimum deblending contrast (default: 0.005; SExtractor default).
    detect_thresh : float
        Detection threshold in sigma (default: 1.0).
    flag : str
        Label appended to auxiliary output filenames.
    gain_key : str | None
        FITS header keyword for detector gain (e⁻/ADU).  None → gain = 1.
    logger : logging.Logger | None
        Logger instance.
    output_dir : str | Path | None
        Directory for the check-image (background-subtracted FITS).
    phot_apertures : np.ndarray | None
        Aperture **diameters** in pixels.  None → single 10-px aperture.
    ref_image : str | Path
        Reference image for dual-image mode (detection on ref, phot on sci).

    Returns
    -------
    astropy.table.Table
        Source table with columns:
        XWIN_IMAGE, YWIN_IMAGE, ALPHAWIN_J2000, DELTAWIN_J2000,
        FLAGS, A_IMAGE, B_IMAGE, THETA_IMAGE, FWHM_IMAGE, FWHM_WORLD,
        KRON_RADIUS, FLUX_RADIUS, FLUX_AUTO, FLUXERR_AUTO,
        MAG_AUTO, MAGERR_AUTO, FLUX_PETRO, FLUXERR_PETRO,
        MAG_PETRO, MAGERR_PETRO,
        MAG_APER[_N], MAGERR_APER[_N], FLUX_APER[_N], FLUXERR_APER[_N],
        [VECTOR_ASSOC, NUMBER_ASSOC when assoc_file is provided].
    """
    fits_path  = Path(fits_path)
    out_dir    = Path(output_dir) if output_dir else fits_path.parent

    # 1. Load images -----------------------------------------------------------
    sci_data, sci_header = _load_fits_data(fits_path)

    detect_data = sci_data
    if ref_image:
        detect_data, _ = _load_fits_data(ref_image)

    # 2. Background estimation -------------------------------------------------
    detect_bkg = sep.Background(detect_data,
                                bw=back_size, bh=back_size,
                                fw=back_filtersize, fh=back_filtersize)
    detect_sub = np.ascontiguousarray((detect_data - detect_bkg).astype(np.float64))

    sci_bkg    = sep.Background(sci_data,
                                bw=back_size, bh=back_size,
                                fw=back_filtersize, fh=back_filtersize)
    sci_sub    = np.ascontiguousarray((sci_data - sci_bkg).astype(np.float64))

    # Save background-subtracted image as check file for diagnostics
    try:
        check_path = out_dir / ('check_' + fits_path.name)
        fits.writeto(str(check_path), sci_sub, sci_header, overwrite=True)
    except Exception as exc:
        if logger:
            logger.warning(f'Could not write check image: {exc}')

    # 3. Source extraction -----------------------------------------------------
    # Global-RMS threshold avoids the (thresh × local_rms²) pitfall that
    # arises when passing an error map to sep.extract.
    thresh = max(float(detect_thresh), float(analysis_thresh)) * detect_bkg.globalrms

    sep.set_extract_pixstack(min(2_000_000, max(500_000, detect_sub.size // 3)))
    sep.set_sub_object_limit(262144)

    try:
        objects = sep.extract(
            detect_sub, thresh,
            deblend_nthresh=int(deblend_nthresh),
            deblend_cont=float(deblend_mincont),
            minarea=5,
            filter_type='matched')
    except Exception as exc:
        msg = f'sep.extract failed on {fits_path}: {exc}'
        print(bcolors.FAIL + msg + bcolors.ENDC)
        if logger:
            logger.error(msg)
        return table.Table()

    if len(objects) == 0:
        msg = f'No sources detected in {fits_path}'
        print(bcolors.WARNING + msg + bcolors.ENDC)
        if logger:
            logger.warning(msg)
        return table.Table()

    # 4. WCS sky coordinates ---------------------------------------------------
    try:
        wcs_obj   = wcs.WCS(sci_header)
        pix_scale = float(np.median(proj_plane_pixel_scales(wcs_obj)) * 3600.0)
        ra, dec   = wcs_obj.wcs_pix2world(objects['x'], objects['y'], 0)
    except Exception as exc:
        if logger:
            logger.warning(f'WCS transformation failed: {exc}')
        ra        = np.full(len(objects), np.nan)
        dec       = np.full(len(objects), np.nan)
        pix_scale = 1.0

    # 5. Gain ------------------------------------------------------------------
    gain_val = get_gain(fits_path, gain_key, logger)

    # 6. Kron (AUTO) photometry ------------------------------------------------
    try:
        kronrad, _ = sep.kron_radius(
            sci_sub,
            objects['x'], objects['y'],
            objects['a'], objects['b'], objects['theta'], 6.0)
        kronrad  = np.maximum(kronrad, 0.0)
        r_kron   = np.maximum(KRON_FACTOR * kronrad, KRON_MIN_R)

        flux_auto, fluxerr_auto, flag_auto = sep.sum_ellipse(
            sci_sub,
            objects['x'], objects['y'],
            objects['a'], objects['b'], objects['theta'],
            r_kron, gain=gain_val, subpix=5)
    except Exception as exc:
        if logger:
            logger.warning(f'Kron photometry failed: {exc}')
        flux_auto    = np.zeros(len(objects))
        fluxerr_auto = np.zeros(len(objects))
        flag_auto    = np.zeros(len(objects), dtype=int)
        kronrad      = np.zeros(len(objects))

    # 7. Half-light radius (→ FWHM) -------------------------------------------
    try:
        flux_radius, _ = sep.flux_radius(
            sci_sub,
            objects['x'], objects['y'],
            6.0 * objects['a'],
            0.5,
            normflux=np.abs(flux_auto),
            subpix=5)
        flux_radius = np.maximum(flux_radius, 0.0)
    except Exception as exc:
        if logger:
            logger.warning(f'flux_radius failed: {exc}')
        flux_radius = np.zeros(len(objects))

    # 8. Petrosian photometry (≈ 2 × half-light radius) -----------------------
    try:
        r_petro = np.maximum(PETRO_FACTOR * flux_radius, PETRO_MIN_R)
        flux_petro, fluxerr_petro, _ = sep.sum_ellipse(
            sci_sub,
            objects['x'], objects['y'],
            objects['a'], objects['b'], objects['theta'],
            r_petro, gain=gain_val, subpix=5)
    except Exception as exc:
        if logger:
            logger.warning(f'Petrosian photometry failed: {exc}')
        flux_petro    = flux_auto.copy()
        fluxerr_petro = fluxerr_auto.copy()

    # 9. Fixed-aperture photometry --------------------------------------------
    if phot_apertures is not None:
        aper_diam = np.asarray(phot_apertures, dtype=float)
    else:
        aper_diam = np.array([10.0])

    aper_radii     = aper_diam / 2.0
    flux_aper_all  = []
    fluxerr_aper_all = []

    for r in aper_radii:
        try:
            fl, fle, _ = sep.sum_circle(
                sci_sub,
                objects['x'], objects['y'],
                float(r), gain=gain_val, subpix=5)
        except Exception:
            fl  = np.zeros(len(objects))
            fle = np.zeros(len(objects))
        flux_aper_all.append(fl)
        fluxerr_aper_all.append(fle)

    # 10. Derived quantities --------------------------------------------------
    fwhm_image = 2.0 * flux_radius          # Gaussian approximation
    fwhm_world = fwhm_image * pix_scale     # arcsec

    def _mag(flux: np.ndarray) -> np.ndarray:
        return np.where(flux > 0,
                        -2.5 * np.log10(np.maximum(flux, 1e-30)),
                        np.nan)

    def _magerr(flux: np.ndarray, ferr: np.ndarray) -> np.ndarray:
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(flux > 0,
                            2.5 / np.log(10.) * np.abs(ferr)
                            / np.maximum(flux, 1e-30),
                            np.nan)

    # 11. Assemble output table -----------------------------------------------
    flags  = objects['flag'].astype(int) | flag_auto.astype(int)

    result = table.Table()
    result['XWIN_IMAGE']     = objects['x'] + 1.0     # 0-indexed → 1-indexed
    result['YWIN_IMAGE']     = objects['y'] + 1.0
    result['ALPHAWIN_J2000'] = ra
    result['DELTAWIN_J2000'] = dec
    result['FLAGS']          = flags
    result['A_IMAGE']        = objects['a']
    result['B_IMAGE']        = objects['b']
    result['THETA_IMAGE']    = np.degrees(objects['theta'])
    result['FWHM_IMAGE']     = fwhm_image
    result['FWHM_WORLD']     = fwhm_world
    result['KRON_RADIUS']    = kronrad
    result['FLUX_RADIUS']    = flux_radius
    result['FLUX_AUTO']      = flux_auto
    result['FLUXERR_AUTO']   = fluxerr_auto
    result['MAG_AUTO']       = _mag(flux_auto)
    result['MAGERR_AUTO']    = _magerr(flux_auto, fluxerr_auto)
    result['FLUX_PETRO']     = flux_petro
    result['FLUXERR_PETRO']  = fluxerr_petro
    result['MAG_PETRO']      = _mag(flux_petro)
    result['MAGERR_PETRO']   = _magerr(flux_petro, fluxerr_petro)

    # First aperture has no numeric suffix; subsequent are _1, _2, …
    for i, (fl, fle) in enumerate(zip(flux_aper_all, fluxerr_aper_all)):
        sfx = '' if i == 0 else f'_{i}'
        result[f'FLUX_APER{sfx}']    = fl
        result[f'FLUXERR_APER{sfx}'] = fle
        result[f'MAG_APER{sfx}']     = _mag(fl)
        result[f'MAGERR_APER{sfx}']  = _magerr(fl, fle)

    # 12. ASSOC positional matching -------------------------------------------
    if assoc_file is not None:
        try:
            assoc_data  = ascii.read(str(assoc_file))
            col_idx     = [int(c.strip()) - 1 for c in assoc_cols.split(',')]
            assoc_x     = np.asarray(assoc_data.columns[col_idx[0]], dtype=float)
            assoc_y     = np.asarray(assoc_data.columns[col_idx[1]], dtype=float)
            src_x       = np.asarray(result['XWIN_IMAGE'])
            src_y       = np.asarray(result['YWIN_IMAGE'])

            matched_src   = []
            matched_assoc = []

            for j, (ax, ay) in enumerate(zip(assoc_x, assoc_y)):
                dists   = np.hypot(src_x - ax, src_y - ay)
                nearest = int(np.argmin(dists))
                if dists[nearest] <= assoc_radius:
                    matched_src.append(nearest)
                    matched_assoc.append(j)

            if not matched_src:
                if logger:
                    logger.warning(
                        f'ASSOC: no match within {assoc_radius:.1f} px '
                        f'in {fits_path}')
                return table.Table()

            result = result[matched_src]

            if len(col_idx) > 2:
                vec = np.asarray(assoc_data.columns[col_idx[2]], dtype=float)
                result['VECTOR_ASSOC'] = vec[matched_assoc]
            else:
                result['VECTOR_ASSOC'] = np.zeros(len(matched_src))

            result['NUMBER_ASSOC'] = (np.asarray(matched_assoc, dtype=float) + 1.0)

        except Exception as exc:
            msg = f'ASSOC matching failed: {exc}'
            print(bcolors.WARNING + msg + bcolors.ENDC)
            if logger:
                logger.warning(msg)
            result['VECTOR_ASSOC'] = np.zeros(len(result))
            result['NUMBER_ASSOC'] = np.zeros(len(result))

    return result


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def postprocess_catalog(
    catalog: table.Table,
    verbose: bool = True,
) -> table.Table:
    """Post-process a source catalog: replace negative fluxes and fix errors.

    Sets all negative-flux entries to NaN (object becomes an upper limit),
    then recomputes asymmetric magnitude errors from flux and flux-error
    columns so the catalog is consistent with forced-photometry upper limits.

    Parameters
    ----------
    catalog : astropy.table.Table
        Source catalog with ``FLUX_*`` / ``FLUXERR_*`` column pairs.
    verbose : bool
        Print informational messages (default: True).

    Returns
    -------
    astropy.table.Table
        Modified catalog (modified in-place, also returned for convenience).
    """
    if verbose:
        print(bcolors.WARNING
              + 'Recalculate magnitude errors  [ Δmag = −2.5 log(F ± ΔF) / F ]'
              + bcolors.ENDC)
        print(bcolors.WARNING
              + 'Objects with F ≤ 0 become 3σ upper limits; MAGERR set to NaN.'
              + bcolors.ENDC)

    flux_keys = [k for k in catalog.colnames
                 if k.startswith('FLUX_') and 'RADIUS' not in k]

    for key in flux_keys:
        err_key = key.replace('FLUX_', 'FLUXERR_')
        if err_key not in catalog.colnames:
            continue
        neg = np.asarray(catalog[key]) <= 0.0
        if np.any(neg):
            catalog[key][neg]     = np.nan
            catalog[err_key][neg] = np.nan

    for err_key in [k for k in catalog.colnames if k.startswith('FLUXERR_')]:
        flux_key = err_key.replace('FLUXERR_', 'FLUX_')
        if flux_key not in catalog.colnames:
            continue

        flux    = np.asarray(catalog[flux_key], dtype=float)
        fluxerr = np.asarray(catalog[err_key],  dtype=float)

        with np.errstate(divide='ignore', invalid='ignore'):
            errp = np.where(flux > 0,
                            -2.5 * np.log10(flux - fluxerr)
                            + 2.5 * np.log10(flux), np.nan)
            errm = np.where(flux > 0,
                            -2.5 * np.log10(flux)
                            + 2.5 * np.log10(flux + fluxerr), np.nan)

        base             = err_key.replace('FLUXERR_', '')
        catalog[f'MAGERR_{base}']  = np.abs(errp + errm) / 2.0
        catalog[f'MAGERRP_{base}'] = errp
        catalog[f'MAGERRM_{base}'] = errm

    return catalog
