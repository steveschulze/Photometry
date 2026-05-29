"""
phot_routines.py — Core photometry algorithms for the ENGRAVE pipeline.

Source detection and photometry use sep (Source Extraction and Photometry),
a C-extension library implementing SExtractor-compatible algorithms without
requiring an external binary.

Dependencies: sep, photutils, astropy, numpy, scipy, pysynphot (HST only)
"""

__version__  = "2026-05-29"
__author__   = "Steve Schulze (steve.schulze@weizmann.ac.il)"

from   astroquery.vizier import Vizier
from   astropy import coordinates as coord
from   astropy import stats, table, time, wcs
from   astropy import units as u
from   astropy.io import ascii, fits
from   astropy.wcs.utils import proj_plane_pixel_scales
import cat_tools
from   cat_tools import catalog_prop
import fits_tools
from   matplotlib import pylab as plt
from   matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from   misc import bcolors
import numpy as np
import numpy.lib.recfunctions as rfn
import os
from   pathlib import Path
from   photutils.aperture import (
    CircularAperture, CircularAnnulus, EllipticalAperture,
    aperture_photometry as phot_aperture_photometry
)
from   plotsettings import *
import sep
import stat_tools
import sys
import warnings
from   astropy.utils.exceptions import AstropyUserWarning
from   scipy import stats as stats_scipy
from   scipy import optimize

warnings.filterwarnings("ignore", category=AstropyUserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Aperture photometry (photutils-based, used for forced photometry)
# ---------------------------------------------------------------------------

def aperture_photometry(IMAGE, POSITIONS, RADII, INNERANNULUS, OUTERANNULUS,
                         RMS, GAIN=1, FA=1, ZEROPOINT=0):
    """Perform circular aperture photometry for one or more apertures.

    Uses photutils for aperture sums and local background subtraction.
    Returns magnitudes in the AB system and flux densities in micro-Jy.

    Parameters
    ----------
    IMAGE : 2-D array
        Background-subtracted science image.
    POSITIONS : tuple or array-like
        (x, y) pixel position(s) of the target (1-indexed, FITS convention).
    RADII : array-like
        Aperture radii in pixels.
    INNERANNULUS, OUTERANNULUS : array-like
        Inner and outer radii of the background annuli in pixels.
    RMS : float
        Global image RMS (used as fallback; local RMS preferred).
    GAIN : float
        Effective detector gain in e⁻/ADU.
    FA : float
        Correlated-noise correction factor (default: 1 = uncorrelated).
    ZEROPOINT : array-like
        AB magnitude zeropoint for each aperture.

    Returns
    -------
    astropy.table.Table
        Photometry table with flux, magnitude, and noise columns per aperture.
    """
    # Pixel positions are 1-indexed (FITS); photutils expects 0-indexed
    pos_0indexed = (np.array(POSITIONS) - 1) if np.ndim(POSITIONS) == 1 else \
                   np.array(POSITIONS) - 1

    src_apers = [CircularAperture(pos_0indexed, r=radius) for radius in RADII]
    bkg_apers = [CircularAnnulus(pos_0indexed, r_in=INNERANNULUS[i],
                                  r_out=OUTERANNULUS[i]) for i in range(len(RADII))]

    src_phot_table = phot_aperture_photometry(IMAGE, src_apers)
    bkg_phot_table = phot_aperture_photometry(IMAGE, bkg_apers)

    # Rename single-aperture column to numbered convention
    if 'aperture_sum_0' not in src_phot_table.colnames:
        src_phot_table.rename_column('aperture_sum', 'aperture_sum_0')
        bkg_phot_table.rename_column('aperture_sum', 'aperture_sum_0')

    # Local background RMS per aperture
    local_rms = [background_local(IMAGE, bkg_apers[i])[-1] for i in range(len(RADII))]

    # Area ratio for background scaling
    area_ratio = [src_apers[i].area / bkg_apers[i].area for i in range(len(RADII))]

    for i in range(len(RADII)):
        bkg_scaled = bkg_phot_table['aperture_sum_' + str(i)] * area_ratio[i]
        src_phot_table['bkg_' + str(i)] = bkg_scaled
        src_phot_table['bkg_sub_aperture_sum_' + str(i)] = \
            src_phot_table['aperture_sum_' + str(i)] - bkg_scaled

    for i in range(len(RADII)):
        src_phot_table['ZP_APER_' + str(i)] = ZEROPOINT[i]

        net_flux = src_phot_table['bkg_sub_aperture_sum_' + str(i)]

        shot_noise = np.sqrt(np.abs(net_flux) / GAIN)
        bkg_noise  = np.sqrt(local_rms[i]**2 * src_apers[i].area / FA)

        src_phot_table['SHOT_NOISE_'   + str(i)] = shot_noise
        src_phot_table['BKG_NOISE_'    + str(i)] = bkg_noise
        src_phot_table['TOTAL_ERROR_'  + str(i)] = np.sqrt(shot_noise**2 + bkg_noise**2)

        # Flux density in micro-Jy (AB zeropoint 23.9 corresponds to 1 µJy)
        factor = 10**(-0.4 * (ZEROPOINT[i] - 23.9))
        fnu     = net_flux * factor
        fnu_err = src_phot_table['TOTAL_ERROR_' + str(i)] * factor

        src_phot_table['FNU_APER_'    + str(i)] = fnu
        src_phot_table['FNUERR_APER_' + str(i)] = fnu_err

        # AB magnitudes
        mag = np.where(fnu > 0,
                       -2.5 * np.log10(fnu) + 23.9,
                       -2.5 * np.log10(np.maximum(3 * fnu_err, 1e-30)) + 23.9)

        mag_errp = np.where(
            fnu > 0,
            -2.5 * np.log10(fnu - fnu_err) + 2.5 * np.log10(fnu),
            np.nan)
        mag_errm = np.where(
            fnu > 0,
            -2.5 * np.log10(fnu) + 2.5 * np.log10(fnu + fnu_err),
            np.nan)

        src_phot_table['MAG_APER_'    + str(i)] = mag
        src_phot_table['MAGERRP_APER_' + str(i)] = mag_errp
        src_phot_table['MAGERRM_APER_' + str(i)] = mag_errm

    for key in src_phot_table.colnames:
        if 'MAG_APER' in key or 'aperture' in key:
            src_phot_table[key].info.format = '%.3f'
        elif 'FNU' in key or 'NOISE' in key or 'TOTAL' in key:
            src_phot_table[key].info.format = '%.3e'

    return src_phot_table


# ---------------------------------------------------------------------------
# Background estimation
# ---------------------------------------------------------------------------

def background(IMAGE, BACK_SIZE=64, BACK_FILTERSIZE=3):
    """Estimate the global RMS of a science image using sep background.

    Parameters
    ----------
    IMAGE : 2-D array
        Science image data.
    BACK_SIZE : int
        Background mesh size in pixels.
    BACK_FILTERSIZE : int
        Background filter size.

    Returns
    -------
    float
        Global background RMS.
    """
    data = np.ascontiguousarray(IMAGE.astype(np.float64))
    bkg  = sep.Background(data, bw=BACK_SIZE, bh=BACK_SIZE,
                          fw=BACK_FILTERSIZE, fh=BACK_FILTERSIZE)
    return bkg.globalrms


def background_local(IMAGE, APERTURE):
    """Compute sigma-clipped statistics inside a background annulus aperture.

    Parameters
    ----------
    IMAGE : 2-D array
        Science image data.
    APERTURE : photutils aperture
        Background annulus aperture object.

    Returns
    -------
    tuple
        (mean, median, std) of pixel values inside the aperture.
    """
    bkg_mask = APERTURE.to_mask(method='center')
    bkg_data = bkg_mask.multiply(IMAGE)
    aper_data = bkg_data[bkg_mask.data > 0]
    return stats.sigma_clipped_stats(aper_data)


# ---------------------------------------------------------------------------
# Gain extraction
# ---------------------------------------------------------------------------

def get_gain(FITS, KEYWORD, LOGGER):
    """Read detector gain from the FITS header.

    Searches for common gain keywords (CCDGAIN, ADCGAIN, ATODGAIN).
    Returns 1 if KEYWORD is None (un-calibrated mode).

    Parameters
    ----------
    FITS : str
        Path to the FITS file.
    KEYWORD : str or None
        Header keyword name, or None to use gain=1.
    LOGGER : logging.Logger
        Logger for error reporting.

    Returns
    -------
    float
        Gain value in e⁻/ADU.
    """
    if KEYWORD is None:
        return 1

    try:
        with fits.open(FITS) as hdu:
            header = hdu[0].header
            if len(hdu) > 1:
                header = hdu[0].header + hdu[1].header
            if len(hdu) > 2:
                header = hdu[0].header + hdu[1].header + hdu[2].header

        for key in ('CCDGAIN', 'ADCGAIN', 'ATODGAIN', KEYWORD):
            if key in header:
                return float(header[key])

        msg = f'Gain keyword ({KEYWORD}) not found in {FITS}'
        print(bcolors.BOLD + bcolors.FAIL + f'\n{msg}\n' + bcolors.ENDC)
        if LOGGER:
            LOGGER.error(msg)
        sys.exit(1)

    except Exception as exc:
        msg = f'Failed to read gain from {FITS}: {exc}'
        print(bcolors.BOLD + bcolors.FAIL + f'\n{msg}\n' + bcolors.ENDC)
        if LOGGER:
            LOGGER.error(msg)
        sys.exit(1)


# ---------------------------------------------------------------------------
# HST photometry
# ---------------------------------------------------------------------------

def hst_aperture_photometry(FITS, POSITIONS, RADII, INNERANNULUS,
                             OUTERANNULUS, PIX2ARCSEC, RMS, FA=1):
    """Aperture photometry for HST drizzled images.

    Computes the effective gain from EXPTIME and CCD gain keyword, applies
    the correlated-noise correction factor for drizzled data, and delegates
    to :func:`aperture_photometry`.

    Parameters
    ----------
    FITS : str
        HST drizzled FITS file path.
    POSITIONS : array-like
        (x, y) target pixel positions (1-indexed).
    RADII : array-like
        Aperture radii in arcseconds.
    INNERANNULUS, OUTERANNULUS : array-like
        Inner/outer background annulus radii in arcseconds.
    PIX2ARCSEC : float
        Pixel scale in arcsec/pixel.
    RMS : float
        Global image RMS.
    FA : float
        Correlated-noise correction (computed internally if not overridden).

    Returns
    -------
    astropy.table.Table
        Photometry table (see :func:`aperture_photometry`).
    """
    try:
        import pysynphot as pyS
    except ImportError as exc:
        raise ImportError(
            "pysynphot is required for HST photometry. "
            "Install with: pip install pysynphot"
        ) from exc

    with fits.open(FITS) as hdulist:
        hdu_header = hdulist[0].header
        if len(hdulist) > 1:
            try:
                hdulist[1].shape
                hdu_header += hdulist[1].header
                hdu_data    = hdulist[1].data
            except Exception:
                hdu_header += hdulist[1].header
                hdu_data    = hdulist[0].data
        else:
            hdu_data = hdulist[0].data

    radii_px          = np.asarray(RADII) / PIX2ARCSEC
    innerannulus_px   = np.asarray(INNERANNULUS) / PIX2ARCSEC
    outerannulus_px   = np.asarray(OUTERANNULUS) / PIX2ARCSEC

    zeropoint_arr     = hst_zeropoint(hdu_header, 2 * np.asarray(RADII))

    # Effective gain: astrodrizzle images are in e⁻/s
    for key in ('CCDGAIN', 'ADCGAIN', 'ATODGAIN'):
        if key in hdu_header:
            gain_key = key
            break
    else:
        raise KeyError(f'No gain keyword found in header of {FITS}')

    effective_gain = hdu_header['EXPTIME'] * hdu_header[gain_key]

    # Correlated-noise correction for drizzled images
    pixfrac = hdu_header['D001PIXF']
    try:
        native_scale = hst_scale(hdu_header['INSTRUME'], hdu_header['APERTURE'])
    except KeyError:
        native_scale = hst_scale(hdu_header['INSTRUME'], None)

    scale = PIX2ARCSEC / native_scale
    fa    = ((scale / pixfrac) * (1. - scale / 3. / pixfrac))**2 \
            if scale < pixfrac else (1. - pixfrac / 3. / scale)**2

    return aperture_photometry(hdu_data, POSITIONS, radii_px,
                                innerannulus_px, outerannulus_px,
                                RMS, GAIN=effective_gain,
                                ZEROPOINT=zeropoint_arr, FA=fa)


def hst_make_cutout(FITS, COORD_OBS, COORD_EXP, RADII, RADII_INNERANNULUS,
                    RADII_OUTERANNULUS, PIX2ARCSEC, OUTDIR):
    """Create a diagnostic multi-panel aperture cutout image for HST data.

    Each panel shows the target region at a different aperture diameter,
    with the aperture and background annulus circles overlaid.

    Parameters
    ----------
    FITS : str
        HST drizzled FITS file path.
    COORD_OBS : tuple
        Observed (x, y) centroid position (1-indexed pixels).
    COORD_EXP : tuple
        Expected (x, y) position from WCS (1-indexed pixels).
    RADII : array-like
        Aperture radii in arcseconds.
    RADII_INNERANNULUS, RADII_OUTERANNULUS : array-like
        Background annulus radii in arcseconds.
    PIX2ARCSEC : float
        Pixel scale in arcsec/pixel.
    OUTDIR : str
        Output directory for the PDF plot.
    """
    with fits.open(FITS) as hdulist:
        header = hdulist[0].header
        if len(hdulist) > 1:
            try:
                hdulist[1].shape
                image = hdulist[1].data
            except Exception:
                image = hdulist[0].data
        else:
            image = hdulist[0].data

    halfwidth   = 100
    cx, cy      = int(COORD_OBS[0]), int(COORD_OBS[1])
    ny, nx      = image.shape
    xmin        = max(cx - halfwidth, 0)
    xmax        = min(cx + halfwidth, nx)
    ymin        = max(cy - halfwidth, 0)
    ymax        = min(cy + halfwidth, ny)

    cutout      = image[ymin:ymax, xmin:xmax]
    flat        = cutout.flatten()
    flat        = flat[np.isfinite(flat) & (flat > 0)]

    vmin        = np.percentile(flat, 55) if len(flat) > 0 else 1e-3
    vmax        = np.percentile(flat, 99) if len(flat) > 0 else 1.0

    num_apers   = len(RADII)
    ncols       = 3
    nrows       = int(np.ceil(num_apers / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(np.sqrt(2) * 9, 9))
    axes_flat = np.array(axes).flatten()

    for i in range(num_apers):
        ax = axes_flat[i]
        ax.imshow(cutout, cmap='binary', interpolation='nearest',
                  origin='lower', norm=LogNorm(vmin=vmin, vmax=vmax))

        ax.plot(halfwidth, halfwidth, marker='x', mew=5, color=color_green, ms=12)

        if len(COORD_OBS) > 0:
            ox = COORD_OBS[0] - cx + halfwidth - 1
            oy = COORD_OBS[1] - cy + halfwidth - 1
            ax.errorbar(ox, oy, mew=5, marker='x', color=vigit_color_1, ms=12)

        for r, color in [(RADII[i] / PIX2ARCSEC, vigit_color_1),
                         (RADII_INNERANNULUS[i] / PIX2ARCSEC, color_yellow),
                         (RADII_OUTERANNULUS[i] / PIX2ARCSEC, color_yellow)]:
            circle = plt.Circle((halfwidth, halfwidth), r,
                                 ls='solid', lw=5, fc='None', ec=color)
            ax.add_patch(circle)

        ax.set_xlim(0, 2 * halfwidth)
        ax.set_ylim(0, 2 * halfwidth)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        for line in ax.yaxis.get_majorticklines() + ax.xaxis.get_majorticklines():
            line.set_markersize(0)

        diam_str = f"\\textbf{{Diameter:}} $\\mathbf{{ {2*RADII[i]:.1f}''}}$"
        ax.text(0.05, 0.95, diam_str, ha='left', va='top',
                transform=ax.transAxes, fontsize=legend_size - 4,
                path_effects=[PathEffects.withStroke(linewidth=6, foreground="w")])

    for j in range(num_apers, len(axes_flat)):
        axes_flat[j].set_visible(False)

    ax_last = axes_flat[num_apers - 1]
    ax_last.text(0.95, 0.175, "\\textbf{Host centroid}", ha='right', va='bottom',
                  transform=ax_last.transAxes, color=color_green,
                  fontsize=legend_size - 4,
                  path_effects=[PathEffects.withStroke(linewidth=6, foreground="w")])
    ax_last.text(0.95, 0.05, "\\textbf{Transient position}", ha='right', va='bottom',
                  transform=ax_last.transAxes, color=vigit_color_1,
                  fontsize=legend_size - 4,
                  path_effects=[PathEffects.withStroke(linewidth=6, foreground="w")])

    plt.savefig(Path(OUTDIR) / (Path(FITS).stem + '.pdf'), dpi=600)
    plt.close()


def hst_cog(FITS, POSITIONS, INNERANNULUS, OUTERANNULUS, PIX2ARCSEC, RMS, OUTDIR):
    """Compute and plot the curve of growth for an HST source.

    Measures the source brightness at 20 logarithmically-spaced aperture
    diameters from 0.2 to 5 arcsec.  Monte Carlo resampling propagates
    asymmetric uncertainties to derive median brightness and scatter.

    Parameters
    ----------
    FITS : str
        HST drizzled FITS file path.
    POSITIONS : array-like
        (x, y) target pixel position (1-indexed).
    INNERANNULUS, OUTERANNULUS : float
        Scaling factors for the background annulus relative to aperture radius.
    PIX2ARCSEC : float
        Pixel scale in arcsec/pixel.
    RMS : float
        Global image RMS.
    OUTDIR : str
        Output directory for plots and data tables.
    """
    apertures    = np.linspace(0.2, 5, 20)
    inner        = INNERANNULUS * apertures
    outer        = OUTERANNULUS * apertures

    photometry   = hst_aperture_photometry(FITS, POSITIONS,
                                            apertures / 2., inner / 2., outer / 2.,
                                            PIX2ARCSEC, RMS)

    mags         = np.array([photometry['MAG_APER_' + str(i)][0]
                              for i in range(len(apertures))])
    mags_errp    = np.array([photometry['MAGERRP_APER_' + str(i)][0]
                              for i in range(len(apertures))])
    mags_errm    = np.array([photometry['MAGERRM_APER_' + str(i)][0]
                              for i in range(len(apertures))])

    mask_det     = np.isfinite(mags_errp) & (mags_errp != 0.)
    mask_ul      = ~mask_det

    # Save raw CoG data
    cog_data = table.Table(
        [apertures, np.round(mags, 3), np.round(mags_errp, 3), np.round(mags_errm, 3)],
        names=('DIAMETER', 'MAG', 'MAGERRP', 'MAGERRM'))
    outdir_path = Path(OUTDIR)
    base        = Path(FITS).stem
    ascii.write(cog_data, outdir_path / (base + '_cog_data.ascii'),
                overwrite=True, format='no_header')

    # Monte Carlo error propagation (vectorized)
    niter        = 10000
    mags_d       = mags[mask_det]
    errp_d       = mags_errp[mask_det]
    errm_d       = mags_errm[mask_det]

    u_mc         = np.random.uniform(size=(niter, len(mags_d)))
    sign         = u_mc < 0.5
    values       = np.where(sign,
                            stats_scipy.norm.ppf(u_mc, mags_d, errm_d),
                            stats_scipy.norm.ppf(u_mc, mags_d, errp_d))
    cog_stats    = stats.sigma_clipped_stats(values, axis=0)[1]  # median per MC iter

    cog_median   = float(np.median(cog_stats))
    cog_std      = float(np.std(cog_stats))

    cog_stats_table = table.Table(
        [[cog_median], [cog_median], [cog_std], [niter]],
        names=('MEAN', 'MEDIAN', 'STD', 'NITER'))
    ascii.write(cog_stats_table, outdir_path / (base + '_cog_stat.ascii'),
                overwrite=True)

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(9 * np.sqrt(2.), 9))

    for sigma, shade in [(3, '0.95'), (2, '0.85'), (1, '0.75')]:
        ax.axhspan(cog_median - sigma * cog_std,
                   cog_median + sigma * cog_std, color=shade)
    ax.axhline(cog_median, lw=3, color='0.75')

    ax.errorbar(apertures, mags, mags_errp, ms=0, lw=3, color=vigit_color_12)
    ax.errorbar(apertures[mask_det], mags[mask_det],
                [mags_errp[mask_det], mags_errm[mask_det]],
                marker='o', ms=9, lw=0, color=vigit_color_12,
                capsize=0, elinewidth=2)
    if mask_ul.any():
        ax.errorbar(apertures[mask_ul], mags[mask_ul],
                    marker='v', ms=9, lw=0, color=vigit_color_12,
                    mec=vigit_color_12)

    ax.text(0.05, 0.95,
            f'median = {cog_median:.3f}, std = {cog_std:.3f}',
            ha='left', va='top', fontsize=label_size, transform=ax.transAxes,
            color='k')
    ax.set_xlabel("Diameter (arcsec)")
    ax.set_ylabel('Brightness (mag, AB)')
    ax.set_xlim(0, apertures[-1])
    valid_mags = mags[np.isfinite(mags)]
    if len(valid_mags):
        ax.set_ylim(np.nanmax(mags), np.nanmin(mags) - 0.5)

    plt.savefig(outdir_path / (base + '_cog.pdf'))
    plt.close()


def hst_scale(INSTRUMENT, MODE):
    """Return the native (pre-drizzle) pixel scale for an HST instrument.

    Parameters
    ----------
    INSTRUMENT : str
        HST instrument name (ACS, WFC3, WFPC2, NICMOS).
    MODE : str or None
        Detector / aperture name (e.g., 'WFC', 'UVIS', 'IR').

    Returns
    -------
    float
        Native pixel scale in arcsec/pixel.
    """
    scales = {
        'NICMOS': 0.05071,
        'WFPC2':  0.10,
    }
    if INSTRUMENT in scales:
        return scales[INSTRUMENT]

    if INSTRUMENT == 'ACS':
        if MODE and 'WFC' in MODE:
            return 0.05
        if MODE == 'HRC':
            return 0.0265
        if MODE == 'SBC':
            return 0.032
        raise ValueError(f'Unknown ACS mode: {MODE}')

    if INSTRUMENT == 'WFC3':
        if MODE and 'UVIS' in MODE:
            return 0.040
        if MODE and 'IR' in MODE:
            return 0.13
        raise ValueError(f'Unknown WFC3 mode: {MODE}')

    raise ValueError(f'Unknown HST instrument: {INSTRUMENT}')


def hst_zeropoint(HEADER, DIAMETER):
    """Compute HST AB zeropoint with pysynphot aperture correction.

    Uses PHOTFLAM and PHOTPLAM from the FITS header to derive the infinite-
    aperture zeropoint, then calls pysynphot to compute the enclosed-flux
    fraction at the requested aperture diameter(s).

    Parameters
    ----------
    HEADER : fits.Header
        FITS header containing HST photometric keywords.
    DIAMETER : array-like
        Aperture diameters in arcseconds. Diameters >= 4'' use no correction.

    Returns
    -------
    np.ndarray
        AB zeropoints (one per aperture diameter) including aperture correction.
    """
    try:
        import pysynphot as pyS
    except ImportError as exc:
        raise ImportError(
            "pysynphot is required for HST zeropoints. "
            "Install with: pip install pysynphot"
        ) from exc

    zp_inf = (-2.5 * np.log10(HEADER['PHOTFLAM'])
              - 5. * np.log10(HEADER['PHOTPLAM'])
              - 2.408)

    spec_bb  = pyS.BlackBody(10000)
    DIAMETER = np.atleast_1d(np.where(np.array(DIAMETER) < 4, DIAMETER, 4.0))

    # Determine filter name
    try:
        filt = HEADER['FILTER1'] if 'CLEAR' in HEADER['FILTER2'] else HEADER['FILTER2']
    except KeyError:
        filt = HEADER.get('FILTER', HEADER.get('FILTNAM1', 'CLEAR'))

    ap_correction = np.ones(len(DIAMETER))

    if HEADER.get('INSTRUME', '').upper() not in ('WFPC2', 'NICMOS'):
        try:
            photmode = HEADER['PHOTMODE'].replace(' ', ',')
            bandprops = [
                pyS.ObsBandpass(f'{photmode},aper#{d:.2f}')
                for d in DIAMETER
            ]
            bandref = pyS.ObsBandpass(f'{photmode},aper#4.00')
        except Exception:
            mjd = int(time.Time(HEADER['DATE-OBS'], format='isot',
                                scale='utc').jd)
            inst = HEADER['INSTRUME']
            det  = HEADER['APERTURE']
            bandprops = [
                pyS.ObsBandpass(f'{inst},{det},{filt},mjd#{mjd},aper#{d:.2f}')
                for d in DIAMETER
            ]
            bandref = pyS.ObsBandpass(
                f'{inst},{det},{filt},mjd#{mjd},aper#4.00')

        count_ref = pyS.Observation(spec_bb, bandref).countrate()
        for i, (d, bp) in enumerate(zip(DIAMETER, bandprops)):
            if d < 4.:
                ap_correction[i] = (pyS.Observation(spec_bb, bp).countrate()
                                    / count_ref)
    else:
        print(bcolors.WARNING
              + f'Aperture correction for {HEADER["INSTRUME"]} not in pysynphot.'
              + '\nThe AP correction is ~−0.1 mag for r=0.5\". Add by hand.'
              + bcolors.ENDC)

    ap_correction_mag = -2.5 * np.log10(ap_correction)
    return zp_inf - ap_correction_mag


# ---------------------------------------------------------------------------
# Local sequence selection
# ---------------------------------------------------------------------------

def local_sequence(CAT, AUTO=False, FILENAME=None, FITS='', LOGGER=None,
                   LOWER=10, PATH='', UPPER=90):
    """Select comparison stars for building the photometric local sequence.

    In AUTO mode the magnitude range is determined from a lookup table of
    telescope/filter combinations.  In interactive mode the user is prompted
    to apply custom magnitude cuts.

    Parameters
    ----------
    CAT : astropy.table.Table
        Cross-matched catalog with columns MAG_INS, MAG_CAT, MAGERR_CAT,
        MAGERR_INS, XWIN_IMAGE, YWIN_IMAGE, ALPHAWIN_J2000, DELTAWIN_J2000.
    AUTO : bool
        Use automatic magnitude range from instrument lookup table.
    FILENAME : str
        Output file path for the cleaned local sequence catalog.
    FITS : str
        Science FITS filename (used to parse telescope/filter for AUTO mode).
    LOGGER : logging.Logger, optional
        Logger instance.
    LOWER : float
        Lower percentile cut on instrumental magnitude (manual mode).
    PATH : str
        Directory for diagnostic plot output.
    UPPER : float
        Upper percentile cut on instrumental magnitude (manual mode).

    Returns
    -------
    dict
        {'NUMSTARS': int, 'CAT': astropy.table.Table}
    """
    if FILENAME is None:
        print(bcolors.FAIL + 'Filename of the local sequence not specified'
              + bcolors.ENDC)
        sys.exit(1)

    instrument_settings = table.Table(
        names=('TELESCOPE', 'FILTER', 'MAG_BRIGHT', 'MAG_FAINT'),
        dtype=('U100', 'U100', 'f', 'f'))
    for row in [
        ('PanSTARRS', 'g', -18, -13), ('PanSTARRS', 'r', -18, -13),
        ('PanSTARRS', 'i', -18, -13), ('PanSTARRS', 'z', -18, -12),
        ('PanSTARRS', 'y', -18, -12),
        ('2MASS', 'J', -10, -6), ('2MASS', 'H', -10, -6),
        ('2MASS', 'K', -10, -6.5),
        ('UKIDSS', 'J', -16, -11), ('UKIDSS', 'H', -16, -11),
        ('UKIDSS', 'K', -16, -11), ('UKIDSS', 'Y', -16, -8),
        ('SDSS', 'u', -7, -3), ('SDSS', 'g', -10, -3),
        ('SDSS', 'r', -10, -3), ('SDSS', 'i', -10, -4),
        ('SDSS', 'z', -10, -4),
    ]:
        instrument_settings.add_row(row)

    def _fit_zp(mag_ins, mag_cat, mag_cat_err):
        """Fit a constant ZP offset between instrumental and catalog mags."""
        pinit   = 0.
        fitfunc = lambda p, x: p + x
        if np.all(mag_cat_err != 0):
            errfunc = lambda p, x, y, e: (y - fitfunc(p, x)) / e
            out = optimize.leastsq(errfunc, pinit,
                                   args=(mag_ins, mag_cat, mag_cat_err),
                                   full_output=True)
        else:
            errfunc = lambda p, x, y: y - fitfunc(p, x)
            out = optimize.leastsq(errfunc, pinit,
                                   args=(mag_ins, mag_cat), full_output=True)
        return out[0][0]

    def _plot_std(cat_all, cat_sel, zp_offset, bright_cut, faint_cut,
                  fignum, show=False):
        """Diagnostic scatter plot: instrumental vs. catalog magnitude."""
        plt.figure(fignum, figsize=(9 * np.sqrt(2.), 9))
        ax = plt.subplot(111)

        xlim = np.array([min(cat_all['MAG_INS']) - 0.2,
                         max(cat_all['MAG_INS']) + 0.2])
        ax.fill_between(xlim, xlim + zp_offset - 0.1, xlim + zp_offset + 0.1,
                        color=vigit_color_12, alpha=0.25)
        ax.plot(xlim, xlim + zp_offset, lw=2, color=vigit_color_12)

        ax.errorbar(cat_all['MAG_INS'], cat_all['MAG_CAT'],
                    cat_all['MAGERR_CAT'], cat_all['MAGERR_INS'],
                    lw=0, ms=10, marker='o', color='0.75',
                    elinewidth=2, capsize=0)
        ax.errorbar(cat_sel['MAG_INS'], cat_sel['MAG_CAT'],
                    cat_sel['MAGERR_CAT'], cat_sel['MAGERR_INS'],
                    lw=0, ms=10, marker='o', color='k',
                    elinewidth=2, capsize=0)

        if bright_cut is not None:
            ax.axvline(bright_cut, color='k', ls='--')
        if faint_cut is not None:
            ax.axvline(faint_cut, color='k', ls='--')

        ax.set_xlabel("Instrumental magnitude (mag)")
        ax.set_ylabel("Apparent magnitude (mag)")
        ax.grid(True)

        if len(cat_sel['MAG_CAT']) > 1:
            ax.set_xlim(min(cat_all['MAG_INS']) - 0.2,
                        min(max(cat_all['MAG_INS']), 0) + 0.2)
            ax.set_ylim(min(cat_sel['MAG_CAT']) - 0.2,
                        max(cat_sel['MAG_CAT']) + 0.2)

        plt.savefig(str(Path(PATH) / (Path(FITS).stem + '_std.pdf')),
                    dpi=600)
        if show:
            plt.show()
        plt.close()

    def _write_loc_seq(cat_out, filename):
        ascii.write(
            np.array([cat_out['XWIN_IMAGE'], cat_out['YWIN_IMAGE'],
                      cat_out['ALPHAWIN_J2000'], cat_out['DELTAWIN_J2000'],
                      cat_out['MAG_CAT'], cat_out['MAGERR_CAT']]).T,
            filename,
            names=['XWIN_IMAGE', 'YWIN_IMAGE',
                   'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'MAG', 'MAG_ERR'],
            overwrite=True)

    # ---- AUTO mode -----------------------------------------------------------
    if AUTO:
        parts    = Path(FITS).stem.split('_')
        tel      = parts[1] if len(parts) > 1 else ''
        filt_str = [p for p in parts[2:] if len(p) == 1]
        filt     = filt_str[0] if filt_str else ''

        mask_row = ((instrument_settings['TELESCOPE'] == tel) &
                    (instrument_settings['FILTER'] == filt))
        if not np.any(mask_row):
            raise ValueError(
                f'No magnitude range defined for telescope={tel}, filter={filt}')

        bright = float(instrument_settings['MAG_BRIGHT'][mask_row][0])
        faint  = float(instrument_settings['MAG_FAINT'][mask_row][0])

        print(bcolors.OKGREEN
              + f'AUTO magnitude range: {bright:.2f} to {faint:.2f}'
              + bcolors.ENDC)
        if LOGGER:
            LOGGER.info(f'Magnitude range: {bright:.2f} - {faint:.2f}')

        mask_good = (bright <= CAT['MAG_INS']) & (CAT['MAG_INS'] <= faint)
        cat_sel   = CAT[mask_good]

        zp = _fit_zp(cat_sel['MAG_INS'], cat_sel['MAG_CAT'],
                     cat_sel['MAGERR_CAT'])
        _plot_std(CAT, cat_sel, zp, bright, faint, fignum=2)

        print(bcolors.OKGREEN
              + f'Number of stars: {mask_good.sum()}'
              + bcolors.ENDC)
        _write_loc_seq(cat_sel, FILENAME)
        return {'NUMSTARS': int(mask_good.sum()), 'CAT': cat_sel}

    # ---- Full range (no interactive cut) ------------------------------------
    if LOWER == 0. and UPPER == 100:
        zp = _fit_zp(CAT['MAG_INS'], CAT['MAG_CAT'], CAT['MAGERR_CAT'])
        _plot_std(CAT, CAT, zp, None, None, fignum=2)
        _write_loc_seq(CAT, FILENAME)
        return {'NUMSTARS': len(CAT), 'CAT': CAT}

    # ---- Manual / percentile cut --------------------------------------------
    if len(CAT['MAG_INS']) > 10:
        bright_cut = np.percentile(CAT['MAG_INS'], LOWER)
        faint_cut  = np.percentile(CAT['MAG_INS'], UPPER)
    else:
        bright_cut = float(CAT['MAG_INS'].min())
        faint_cut  = float(CAT['MAG_INS'].max())

    mask_good = ((CAT['MAG_INS'] >= bright_cut) & (CAT['MAG_INS'] <= faint_cut))
    cat_sel   = CAT[mask_good]
    zp        = _fit_zp(cat_sel['MAG_INS'], cat_sel['MAG_CAT'],
                        cat_sel['MAGERR_CAT'])
    _plot_std(CAT, cat_sel, zp, bright_cut, faint_cut, fignum=2, show=True)

    print('Current magnitude cuts:')
    print(f'  Bright: {bright_cut:.2f}')
    print(f'  Faint:  {faint_cut:.2f}')

    var1 = input(bcolors.BOLD + bcolors.OKBLUE
                 + '\nWould you like to apply a magnitude cut? [y|[n]] '
                 + bcolors.ENDC)

    cat_cleaned = []

    if var1 != 'y':
        _write_loc_seq(cat_sel, FILENAME)
    else:
        while var1 == 'y':
            bright_cut = float(input(bcolors.WARNING
                                     + 'Limit for the bright stars: '
                                     + bcolors.ENDC))
            faint_cut  = float(input(bcolors.WARNING
                                     + 'Limit for the faint stars: '
                                     + bcolors.ENDC))

            cat_cleaned = CAT[((CAT['MAG_INS'] >= bright_cut)
                               & (CAT['MAG_INS'] <= faint_cut))]
            print(bcolors.OKGREEN
                  + f'\nNumber of stars: {len(cat_cleaned)}\n'
                  + bcolors.ENDC)

            if len(cat_cleaned):
                zp = _fit_zp(cat_cleaned['MAG_INS'], cat_cleaned['MAG_CAT'],
                             cat_cleaned['MAGERR_CAT'])
            _plot_std(CAT, cat_cleaned, zp, bright_cut, faint_cut,
                      fignum=3, show=True)
            _write_loc_seq(cat_cleaned, FILENAME)

            var1 = input(bcolors.BOLD + bcolors.WARNING
                         + 'Would you like to change the cuts? [y|[n]] '
                         + bcolors.ENDC)

    if LOGGER:
        LOGGER.info(f'Magnitude range: {bright_cut:.2f} - {faint_cut:.2f}')

    final_cat = cat_cleaned if len(cat_cleaned) else cat_sel
    if len(final_cat) == 0:
        print(bcolors.FAIL
              + 'No stars in the cleaned catalog. Check your input parameters.'
              + bcolors.ENDC)
        sys.exit(1)

    return {'NUMSTARS': len(final_cat), 'CAT': final_cat}


# ---------------------------------------------------------------------------
# Poststamp (cutout) visualisation
# ---------------------------------------------------------------------------

def make_poststamp(FITS, COORD_EXP, COORD_OBS, PATH=''):
    """Create a two-panel diagnostic poststamp showing expected vs observed.

    Opens the SExtractor / sep check image (``check_<FITS>``) and the
    science frame.  Marks the expected position (from WCS) and the observed
    centroid.

    Parameters
    ----------
    FITS : str
        Science FITS filename.
    COORD_EXP : tuple
        Expected (x, y) position from WCS (1-indexed pixels).
    COORD_OBS : tuple
        Observed (x, y) centroid from source detection (1-indexed pixels).
    PATH : str
        Directory for the output PDF.
    """
    fits_base = Path(FITS).name
    check_path = Path('check_' + fits_base)

    def _load_image(filepath):
        with fits.open(filepath) as hdu:
            if len(hdu) > 1:
                img    = hdu[1].data
                hdr    = hdu[1].header
            else:
                img    = hdu[0].data
                hdr    = hdu[0].header
        return img, hdr

    def _vrange(img, cx, cy, halfwidth=50):
        ny, nx = img.shape
        xmin   = max(int(cx) - halfwidth, 0)
        xmax   = min(int(cx) + halfwidth, nx)
        ymin   = max(int(cy) - halfwidth, 0)
        ymax   = min(int(cy) + halfwidth, ny)
        cutout = img[ymin:ymax, xmin:xmax]
        flat   = np.array(sorted(cutout.flatten()[~np.isnan(cutout.flatten())]))
        if len(flat) == 0:
            return 1e-3, 1.0
        pcts   = np.array([np.percentile(flat, x) for x in range(50, 99)])
        vmin   = pcts[pcts > 0][1] if np.any(pcts > 0) else flat.max() / 0.5
        vmax   = np.percentile(flat, 99)
        return vmin, vmax

    halfwidth = 50
    cx_exp, cy_exp = int(COORD_EXP[0]), int(COORD_EXP[1])

    if check_path.exists():
        check_img, check_hdr = _load_image(check_path)
    else:
        check_img, check_hdr = _load_image(FITS)

    sci_img, sci_hdr = _load_image(FITS)

    def _slice(img, cx, cy):
        ny, nx = img.shape
        return img[max(cy-halfwidth, 0):min(cy+halfwidth, ny),
                   max(cx-halfwidth, 0):min(cx+halfwidth, nx)]

    check_cut = _slice(check_img, cx_exp, cy_exp)
    sci_cut   = _slice(sci_img,   cx_exp, cy_exp)

    check_vmin, check_vmax = _vrange(check_img, cx_exp, cy_exp, halfwidth)
    sci_vmin,   sci_vmax   = _vrange(sci_img,   cx_exp, cy_exp, halfwidth)

    fig, (ax_check, ax_sci) = plt.subplots(
        1, 2, figsize=(np.sqrt(2) * 9, 9))

    for ax, img, vmin, vmax in [
        (ax_check, check_cut, check_vmin, check_vmax),
        (ax_sci,   sci_cut,   sci_vmin,   sci_vmax)
    ]:
        ax.imshow(img, cmap='binary', interpolation='nearest', origin='lower',
                  norm=LogNorm(vmin=max(vmin, 1e-10), vmax=max(vmax, 1e-9)))
        ax.plot(halfwidth, halfwidth, marker='x', mew=5,
                color=vigit_color_1, ms=12)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        for line in ax.yaxis.get_majorticklines() + ax.xaxis.get_majorticklines():
            line.set_markersize(0)
        ax.set_xlim(0, 2 * halfwidth)
        ax.set_ylim(0, 2 * halfwidth)

    if len(COORD_OBS) > 0 and COORD_OBS[0] is not None:
        ox = COORD_OBS[0] - cx_exp + halfwidth - 1
        oy = COORD_OBS[1] - cy_exp + halfwidth - 1
        for ax in (ax_check, ax_sci):
            ax.errorbar(ox, oy, mew=5, marker='x', color=color_green, ms=12)

    ax_sci.text(0.95, 0.95, "\\textbf{Observed}", ha='right', va='top',
                transform=ax_sci.transAxes, color=color_green,
                fontsize=legend_size,
                path_effects=[PathEffects.withStroke(linewidth=6,
                                                      foreground="w")])
    ax_sci.text(0.95, 0.85, "\\textbf{Expected}", ha='right', va='top',
                transform=ax_sci.transAxes, color=vigit_color_1,
                fontsize=legend_size,
                path_effects=[PathEffects.withStroke(linewidth=6,
                                                      foreground="w")])

    outname = str(Path(PATH) / (Path(FITS).stem + '_poststamp.pdf'))
    plt.savefig(outname, dpi=600)
    plt.close()


# ---------------------------------------------------------------------------
# Output catalog builder
# ---------------------------------------------------------------------------

def make_scicat(FITS, OBJECT_PROPERTIES, SEXTRACTOR_PHOTOMETRY,
                PHOTUTILS_PHOTOMETRY, ZEROPOINT_SUMMARY, OFFSET, LOGGER):
    """Compile the final science photometry output table.

    Merges header metadata, WCS positions, detection status, and calibrated
    magnitudes from all apertures into a single summary table.

    Parameters
    ----------
    FITS : str
        Science FITS filename.
    OBJECT_PROPERTIES : astropy.table.Table
        Table with target RA, DEC, X_EXP, Y_EXP, OBJECT, PHOTCAL.
    SEXTRACTOR_PHOTOMETRY : astropy.table.Table or list
        Detections near the science object position (empty list if none).
    PHOTUTILS_PHOTOMETRY : astropy.table.Table
        Forced photometry results (used when SEXTRACTOR_PHOTOMETRY is empty).
    ZEROPOINT_SUMMARY : astropy.table.Table
        Zeropoint summary table from :func:`zeropoint`.
    OFFSET : float
        Detection search radius in arcseconds.
    LOGGER : logging.Logger
        Logger instance.

    Returns
    -------
    astropy.table.Table
        Summary photometry table with PROPERTY / VALUE / ERROR+/- / COMMENT.
    """
    hdu_header = fits.getheader(FITS)
    catalog    = table.Table(
        names=('PROPERTY', 'VALUE', 'ERROR+', 'ERROR-', 'COMMENT'),
        dtype=('U100', 'f', 'f', 'f', 'U100'))

    catalog.add_row(['FILENAME', np.nan, np.nan, np.nan, FITS])

    for key in ('DATE-OBS', 'MJD', 'EXPTIME', 'NCOMBINE'):
        if key in hdu_header:
            catalog.add_row([key, np.nan, np.nan, np.nan, str(hdu_header[key])])
        else:
            catalog.add_row([key, np.nan, np.nan, np.nan,
                             '1' if key == 'NCOMBINE' else '...'])

    # Handle non-standard date/time keywords
    if 'OBSMJD' in hdu_header:
        catalog['COMMENT'][catalog['PROPERTY'] == 'DATE-OBS'] = \
            time.Time(hdu_header['OBSMJD'], format='mjd').isot
        catalog['COMMENT'][catalog['PROPERTY'] == 'MJD'] = \
            str(hdu_header['OBSMJD'])

    if 'T' not in str(catalog['COMMENT'][catalog['PROPERTY'] == 'DATE-OBS'][0]):
        for date_key, time_key in (('ut-date', 'ut-time'), ('UT-DATE', 'UT-TIME')):
            if date_key in hdu_header and time_key in hdu_header:
                try:
                    utc = f"{hdu_header[date_key]}T{hdu_header[time_key]}"
                    catalog['COMMENT'][catalog['PROPERTY'] == 'DATE-OBS'] = utc
                    mjd = time.Time(utc, format='isot', scale='utc').mjd
                    catalog['COMMENT'][catalog['PROPERTY'] == 'MJD'] = \
                        f'{mjd:.7f}'
                    break
                except Exception:
                    pass

    catalog.add_row(['PHOTCAL', np.nan, np.nan, np.nan,
                     OBJECT_PROPERTIES['PHOTCAL'][0]])
    catalog.add_row(['RA',  float(f"{OBJECT_PROPERTIES['RA'][0]:.6f}"),
                     np.nan, np.nan, 'degree'])
    catalog.add_row(['DEC', float(f"{OBJECT_PROPERTIES['DEC'][0]:.6f}"),
                     np.nan, np.nan, 'degree'])
    catalog.add_row(['X_IMAGE_EXP',
                     float(f"{OBJECT_PROPERTIES['X_EXP'][0]:.1f}"),
                     np.nan, np.nan, 'pixel'])
    catalog.add_row(['Y_IMAGE_EXP',
                     float(f"{OBJECT_PROPERTIES['Y_EXP'][0]:.1f}"),
                     np.nan, np.nan, 'pixel'])

    pix_scale = fits_tools.pix2arcsec(FITS)

    if len(SEXTRACTOR_PHOTOMETRY) > 0:
        msg = (f'One or more objects found within {OFFSET:.1f} arcsec '
               f'from the source position')
        print(bcolors.OKGREEN + f'\n{msg}\n' + bcolors.ENDC)
        if LOGGER:
            LOGGER.info(msg)

        x_exp = float(OBJECT_PROPERTIES['X_EXP'][0])
        y_exp = float(OBJECT_PROPERTIES['Y_EXP'][0])
        distance = [float(np.sqrt((float(x['XWIN_IMAGE']) - x_exp)**2
                                   + (float(x['YWIN_IMAGE']) - y_exp)**2))
                    for x in SEXTRACTOR_PHOTOMETRY]
        SEXTRACTOR_PHOTOMETRY['DISTANCE (px)']     = distance
        SEXTRACTOR_PHOTOMETRY['DISTANCE (arcsec)'] = \
            SEXTRACTOR_PHOTOMETRY['DISTANCE (px)'] * pix_scale
        SEXTRACTOR_PHOTOMETRY.sort('DISTANCE (px)')

        for key in ('XWIN_IMAGE', 'YWIN_IMAGE'):
            catalog.add_row([key + '_OBS',
                             float(f'{SEXTRACTOR_PHOTOMETRY[key][0]:.1f}'),
                             np.nan, np.nan, 'pixel'])
        for key in ('ALPHAWIN_J2000', 'DELTAWIN_J2000'):
            catalog.add_row([key + '_OBS',
                             float(f'{SEXTRACTOR_PHOTOMETRY[key][0]:.5f}'),
                             np.nan, np.nan, 'degree'])

        catalog.add_row(['DISTANCE (px)',
                         float(f'{SEXTRACTOR_PHOTOMETRY["DISTANCE (px)"][0]:.2f}'),
                         np.nan, np.nan, 'pixel'])
        catalog.add_row(['DISTANCE (arcsec)',
                         float(f'{SEXTRACTOR_PHOTOMETRY["DISTANCE (arcsec)"][0]:.2f}'),
                         np.nan, np.nan, 'arcsec'])
        catalog.add_row(['FLUX_RADIUS (arcsec)',
                         float(f'{SEXTRACTOR_PHOTOMETRY["FLUX_RADIUS"][0] * pix_scale:.2f}'),
                         np.nan, np.nan, 'arcsec'])

        for key in [k for k in SEXTRACTOR_PHOTOMETRY.colnames if 'MAG_' in k]:
            errp_key = key.replace('MAG_', 'MAGERRP_')
            errm_key = key.replace('MAG_', 'MAGERRM_')
            catalog.add_row([key,
                             float(f'{SEXTRACTOR_PHOTOMETRY[key][0]:.3f}'),
                             float(f'{SEXTRACTOR_PHOTOMETRY[errp_key][0]:.3f}'),
                             float(f'{SEXTRACTOR_PHOTOMETRY[errm_key][0]:.3f}'),
                             'mag'])
    else:
        msg = (f'No object found within {OFFSET:.1f} arcsec '
               f'from the source position')
        print(bcolors.WARNING + f'\n{msg}\n' + bcolors.ENDC)
        if LOGGER:
            LOGGER.info(msg)

        for key in ('XWIN_IMAGE', 'YWIN_IMAGE'):
            catalog.add_row([key + '_OBS', np.nan, np.nan, np.nan, 'pixel'])
        catalog.add_row(['DISTANCE (px)',     np.nan, np.nan, np.nan, 'pixel'])
        catalog.add_row(['DISTANCE (arcsec)', np.nan, np.nan, np.nan, 'arcsec'])

        mag_keys = [k for k in PHOTUTILS_PHOTOMETRY.colnames
                    if k.startswith('MAG_APER_') and not k.startswith('MAG_APER_PHOTUTILS')]
        # Sort by aperture index
        mag_keys = sorted(mag_keys,
                          key=lambda k: int(k.split('_')[-1])
                          if k.split('_')[-1].isdigit() else 0)

        for i, key in enumerate(mag_keys):
            idx     = key.split('_')[-1]
            errp_key = 'MAGERRP_APER_' + idx
            errm_key = 'MAGERRM_APER_' + idx
            fnu_key  = 'FNU_APER_' + idx
            fnuerr   = float(PHOTUTILS_PHOTOMETRY['FNUERR_APER_' + idx][0])

            suffix   = '' if i == 0 else str(i)
            mag_3sig = 23.9 - 2.5 * np.log10(3 * fnuerr)

            catalog.add_row([f'MAG_APER{suffix}_PHOTUTILS',
                             float(f'{PHOTUTILS_PHOTOMETRY[key][0]:.3f}'),
                             float(f'{PHOTUTILS_PHOTOMETRY[errp_key][0]:.3f}'),
                             float(f'{PHOTUTILS_PHOTOMETRY[errm_key][0]:.3f}'),
                             'mag'])
            catalog.add_row([f'MAG_APER{suffix}_PHOTUTILS_3SIGMA',
                             float(f'{mag_3sig:.3f}'), np.nan, np.nan, 'mag'])
            catalog.add_row([f'FNU_APER{suffix}_PHOTUTILS',
                             float(f'{PHOTUTILS_PHOTOMETRY[fnu_key][0]:.3e}'),
                             fnuerr, fnuerr, 'microJy'])

    return catalog


# ---------------------------------------------------------------------------
# SExtractor config file writer (kept for astrometry.net compatibility)
# ---------------------------------------------------------------------------

def setup_sextractor(outdir='.'):
    """Write default SExtractor configuration files.

    Writes ``default.conv``, ``default.nnw``, ``default.param``, and
    ``default.sex`` to *outdir*.  These are needed by astrometry.net's
    ``solve-field`` command which internally calls SExtractor for source
    detection.

    Parameters
    ----------
    outdir : str
        Directory where the configuration files are written (default: '.').
    """
    outdir = Path(outdir)

    default_conv = """CONV NORM
# 3x3 all-ground convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""
    default_nnw = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
 3 10 10  1
-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00
-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00
-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00
 0.00000e+00
 1.00000e+00
"""
    default_param = """NUMBER
XWIN_IMAGE
YWIN_IMAGE
ALPHAWIN_J2000
DELTAWIN_J2000
FLAGS
FWHM_IMAGE
CLASS_STAR
FLUX_AUTO
FLUXERR_AUTO
FLUX_RADIUS
A_IMAGE
B_IMAGE
THETA_IMAGE
KRON_RADIUS
"""
    default_sex = """CATALOG_NAME     test.cat
CATALOG_TYPE     ASCII_HEAD
PARAMETERS_NAME  default.param
DETECT_TYPE      CCD
DETECT_MINAREA   3
DETECT_THRESH    1.5
ANALYSIS_THRESH  1.5
FILTER           Y
FILTER_NAME      default.conv
DEBLEND_NTHRESH  32
DEBLEND_MINCONT  0.005
CLEAN            Y
CLEAN_PARAM      1.0
WEIGHT_TYPE      NONE
PHOT_APERTURES   5
PHOT_AUTOPARAMS  2.5, 3.5
PHOT_PETROPARAMS 2.0, 3.5
PHOT_FLUXFRAC    0.5
SATUR_LEVEL      50000.0
SATUR_KEY        SATURATE
MAG_ZEROPOINT    0.0
GAIN             0.0
GAIN_KEY         GAIN
PIXEL_SCALE      1.0
SEEING_FWHM      1.2
STARNNW_NAME     default.nnw
BACK_TYPE        AUTO
BACK_SIZE        64
BACK_FILTERSIZE  3
CHECKIMAGE_TYPE  NONE
MEMORY_OBJSTACK  3000
MEMORY_PIXSTACK  300000
MEMORY_BUFSIZE   1024
VERBOSE_TYPE     QUIET
WRITE_XML        N
"""
    for fname, content in (
        ('default.conv',  default_conv),
        ('default.nnw',   default_nnw),
        ('default.param', default_param),
        ('default.sex',   default_sex),
    ):
        (outdir / fname).write_text(content)


# ---------------------------------------------------------------------------
# sep-based source extraction and photometry
# ---------------------------------------------------------------------------

def _load_fits_data(fits_path):
    """Load image data and combined header from a FITS file.

    Handles both single-extension and multi-extension FITS (e.g., HST DRZ).

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.

    Returns
    -------
    tuple
        (data : np.ndarray float64, header : fits.Header)
    """
    with fits.open(fits_path) as hdulist:
        header = hdulist[0].header.copy()
        if len(hdulist) > 1:
            try:
                # Check if extension 1 has image data
                _ = hdulist[1].data.shape
                header += hdulist[1].header
                data    = hdulist[1].data
            except (AttributeError, TypeError):
                data = hdulist[0].data
        else:
            data = hdulist[0].data

    # sep requires native byte order and float64
    data = np.array(data, dtype=np.float64)
    if not data.flags['C_CONTIGUOUS']:
        data = np.ascontiguousarray(data)
    return data, header


def sextractor_photometry(
        ANALYSIS_THRESH=1.0,
        ASSOC_NAME=None,
        ASSOC_PARAMS="1,2",
        ASSOC_RADIUS=10.0,
        BACK_SIZE=64,
        BACK_FILTERSIZE=3,
        DEBLEND_NTHRESH=32,
        DEBLEND_MINCONT=0.005,
        DETECT_THRESH=1.0,
        FITS="",
        FLAG="",
        GAIN=None,
        LOGGER=None,
        LOGLEVEL='WARNING',
        PARAMS="DEFAULT",
        PATH=None,
        PHOT_APERTURES=None,
        REF_FILE=""):
    """Detect sources and measure photometry using sep (SExtractor in Python).

    Replaces the original sewpy/SExtractor backend with the pure-Python sep
    library.  The returned table has the same column structure as the original
    SExtractor output so that all downstream code is unaffected.

    **Coordinate convention:** sep positions are 0-indexed.  All pixel
    columns (XWIN_IMAGE, YWIN_IMAGE) are converted to 1-indexed (FITS/
    SExtractor) convention before return.

    **ASSOC mode:** When *ASSOC_NAME* is provided, only sources within
    *ASSOC_RADIUS* pixels of any input position are returned (same semantics
    as SExtractor's ASSOCSELEC_TYPE=MATCHED).

    **Dual-image mode:** When *REF_FILE* is provided, source positions are
    detected on the reference image and photometry is measured on the science
    image (same as running ``sex REF_FILE,FITS``).

    Parameters
    ----------
    ANALYSIS_THRESH : float
        Analysis threshold in sigma for background map (default: 1.0).
    ASSOC_NAME : str, optional
        Path to ASCII file with x, y positions for positional matching.
        Columns are 1-indexed pixel positions.
    ASSOC_PARAMS : str
        Column indices (1-based) for x, y [, mag] in ASSOC_NAME.
    ASSOC_RADIUS : float
        Matching radius in pixels (default: 10.0).
    BACK_SIZE : int
        Background mesh size (default: 64 pixels).
    BACK_FILTERSIZE : int
        Background filter size (default: 3).
    DEBLEND_NTHRESH : int
        Number of deblending sub-thresholds (default: 64).
    DEBLEND_MINCONT : float
        Minimum contrast for deblending (default: 1e-5).
    DETECT_THRESH : float
        Detection threshold in sigma (default: 1.0).
    FITS : str
        Path to the science FITS file.
    FLAG : str
        Label for auxiliary output files.
    GAIN : str, optional
        Header keyword for detector gain (e⁻/ADU). None → gain=1.
    LOGGER : logging.Logger, optional
        Logger instance for messages.
    LOGLEVEL : str
        Logging level string (default: 'WARNING').
    PARAMS : str
        'DEFAULT' or list of column names (currently only DEFAULT supported).
    PATH : str, optional
        Directory for auxiliary files (check image, catalog).
    PHOT_APERTURES : np.ndarray, optional
        Aperture **diameters** in pixels. If None, a single 10-px aperture
        is used.  Multiple apertures produce indexed column names:
        MAG_APER (first), MAG_APER_1 (second), etc.
    REF_FILE : str
        Reference image for dual-image mode. Empty string = single-image.

    Returns
    -------
    astropy.table.Table
        Source table with columns:
        XWIN_IMAGE, YWIN_IMAGE, ALPHAWIN_J2000, DELTAWIN_J2000,
        FLAGS, A_IMAGE, B_IMAGE, THETA_IMAGE, FWHM_IMAGE, FWHM_WORLD,
        KRON_RADIUS, FLUX_RADIUS, FLUX_AUTO, FLUXERR_AUTO, MAG_AUTO,
        MAGERR_AUTO, FLUX_PETRO, FLUXERR_PETRO, MAG_PETRO, MAGERRR_PETRO,
        MAG_APER[_N], MAGERR_APER[_N], FLUX_APER[_N], FLUXERR_APER[_N],
        [VECTOR_ASSOC, NUMBER_ASSOC if ASSOC_NAME given]
    """
    fits_path = Path(FITS)
    path_out  = Path(PATH) if PATH else fits_path.parent

    # 1. Load image data -------------------------------------------------------
    sci_data, sci_header = _load_fits_data(FITS)

    detect_data = sci_data
    if REF_FILE:
        detect_data, _ = _load_fits_data(REF_FILE)

    # 2. Background estimation -------------------------------------------------
    detect_bkg = sep.Background(detect_data,
                                bw=BACK_SIZE, bh=BACK_SIZE,
                                fw=BACK_FILTERSIZE, fh=BACK_FILTERSIZE)
    detect_sub = (detect_data - detect_bkg).astype(np.float64)
    detect_sub = np.ascontiguousarray(detect_sub)

    sci_bkg    = sep.Background(sci_data,
                                bw=BACK_SIZE, bh=BACK_SIZE,
                                fw=BACK_FILTERSIZE, fh=BACK_FILTERSIZE)
    sci_sub    = (sci_data - sci_bkg).astype(np.float64)
    sci_sub    = np.ascontiguousarray(sci_sub)

    # Save background-subtracted image as check image (replaces SExtractor
    # CHECKIMAGE_TYPE APERTURES for the poststamp diagnostic)
    check_path = path_out / ('check_' + fits_path.name)
    try:
        fits.writeto(str(check_path), sci_sub, sci_header, overwrite=True)
    except Exception as exc:
        if LOGGER:
            LOGGER.warning(f'Could not write check image: {exc}')

    # 3. Source extraction -----------------------------------------------------
    # Use global-RMS-based absolute threshold (consistent with SExtractor's
    # default BACK_TYPE=AUTO without a weight map).
    # Do NOT pass err= here: sep interprets thresh as thresh*err[i,j] when
    # an error map is provided, which gives thresh = N_sigma × local_rms²
    # (far too low) instead of the intended N_sigma × globalrms.
    thresh = max(float(DETECT_THRESH), float(ANALYSIS_THRESH)) * detect_bkg.globalrms

    # Scale sep internal buffers to image size; defaults of 300k pixels /
    # 1024 sub-objects are too small for wide-field frames.
    sep.set_extract_pixstack(min(2_000_000, max(500_000, detect_sub.size // 3)))
    sep.set_sub_object_limit(262144)

    try:
        objects = sep.extract(
            detect_sub, thresh,
            deblend_nthresh=int(DEBLEND_NTHRESH),
            deblend_cont=float(DEBLEND_MINCONT),
            minarea=5,
            filter_type='matched')
    except Exception as exc:
        msg = f'sep.extract failed on {FITS}: {exc}'
        print(bcolors.FAIL + msg + bcolors.ENDC)
        if LOGGER:
            LOGGER.error(msg)
        return table.Table()

    if len(objects) == 0:
        msg = f'No sources detected in {FITS}'
        print(bcolors.WARNING + msg + bcolors.ENDC)
        if LOGGER:
            LOGGER.warning(msg)
        return table.Table()

    # 4. WCS transformation for sky coordinates --------------------------------
    try:
        wcs_obj    = wcs.WCS(sci_header)
        pix_scale  = float(np.median(proj_plane_pixel_scales(wcs_obj)) * 3600.)
        # sep: 0-indexed; WCS origin=0 → 0-indexed input
        ra, dec    = wcs_obj.wcs_pix2world(objects['x'], objects['y'], 0)
    except Exception as exc:
        if LOGGER:
            LOGGER.warning(f'WCS transformation failed: {exc}')
        ra  = np.full(len(objects), np.nan)
        dec = np.full(len(objects), np.nan)
        pix_scale = 1.0

    # 5. Gain ------------------------------------------------------------------
    gain_val = float(get_gain(FITS, GAIN, LOGGER))

    # 6. Kron (AUTO) photometry ------------------------------------------------
    # SExtractor PHOT_AUTOPARAMS = 2.5, 3.5
    kron_factor = 2.5
    kron_min_r  = 3.5

    try:
        kronrad, krflag = sep.kron_radius(
            sci_sub,
            objects['x'], objects['y'],
            objects['a'], objects['b'], objects['theta'],
            6.0)
        kronrad = np.maximum(kronrad, 0.0)
        r_kron  = np.maximum(kron_factor * kronrad, kron_min_r)

        flux_auto, fluxerr_auto, flag_auto = sep.sum_ellipse(
            sci_sub,
            objects['x'], objects['y'],
            objects['a'], objects['b'], objects['theta'],
            r_kron,
            gain=gain_val, subpix=5)
    except Exception as exc:
        if LOGGER:
            LOGGER.warning(f'Kron photometry failed: {exc}')
        flux_auto    = np.zeros(len(objects))
        fluxerr_auto = np.zeros(len(objects))
        flag_auto    = np.zeros(len(objects), dtype=int)
        kronrad      = np.zeros(len(objects))

    # 7. Half-light radius (FLUX_RADIUS at 50% enclosed flux) -----------------
    try:
        # Reference aperture: 6 * semi-major axis
        flux_radius, fr_flag = sep.flux_radius(
            sci_sub,
            objects['x'], objects['y'],
            6. * objects['a'],
            0.5,
            normflux=np.abs(flux_auto),
            subpix=5)
        flux_radius = np.maximum(flux_radius, 0.0)
    except Exception as exc:
        if LOGGER:
            LOGGER.warning(f'flux_radius computation failed: {exc}')
        flux_radius = np.zeros(len(objects))

    # 8. Petrosian photometry (approximation: 2.0 × half-light radius) --------
    # SExtractor PHOT_PETROPARAMS = 2.0, 3.5
    petro_factor = 2.0
    petro_min_r  = 3.5

    try:
        r_petro = np.maximum(petro_factor * flux_radius, petro_min_r)
        flux_petro, fluxerr_petro, flag_petro = sep.sum_ellipse(
            sci_sub,
            objects['x'], objects['y'],
            objects['a'], objects['b'], objects['theta'],
            r_petro,
            gain=gain_val, subpix=5)
    except Exception as exc:
        if LOGGER:
            LOGGER.warning(f'Petrosian photometry failed: {exc}')
        flux_petro    = flux_auto.copy()
        fluxerr_petro = fluxerr_auto.copy()

    # 9. Fixed-aperture photometry --------------------------------------------
    if isinstance(PHOT_APERTURES, (np.ndarray, list)):
        aper_diameters = np.asarray(PHOT_APERTURES, dtype=float)
    else:
        aper_diameters = np.array([10.0])  # 10-pixel diameter default

    aper_radii = aper_diameters / 2.0

    flux_aper_all    = []
    fluxerr_aper_all = []
    for r in aper_radii:
        try:
            fl, fle, _ = sep.sum_circle(
                sci_sub,
                objects['x'], objects['y'],
                r, gain=gain_val, subpix=5)
        except Exception:
            fl  = np.zeros(len(objects))
            fle = np.zeros(len(objects))
        flux_aper_all.append(fl)
        fluxerr_aper_all.append(fle)

    # 10. Derived quantities --------------------------------------------------
    # FWHM from half-light radius (Gaussian approximation):
    # For a Gaussian, r_half ≈ 1.177 σ and FWHM = 2.355 σ → FWHM ≈ 2 r_half
    fwhm_image = 2. * flux_radius
    fwhm_world = fwhm_image * pix_scale

    def _safe_mag(flux):
        """Convert flux to instrumental magnitude (ZP = 0)."""
        return np.where(flux > 0, -2.5 * np.log10(np.maximum(flux, 1e-30)), np.nan)

    def _safe_magerr(flux, fluxerr):
        """Flux-ratio magnitude error (symmetric approximation)."""
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(flux > 0,
                            2.5 / np.log(10.) * np.abs(fluxerr) / np.maximum(flux, 1e-30),
                            np.nan)

    # 11. Assemble output table -----------------------------------------------
    flags = objects['flag'].astype(int) | flag_auto.astype(int)

    result = table.Table()
    result['XWIN_IMAGE']     = objects['x'] + 1.  # 0→1-indexed
    result['YWIN_IMAGE']     = objects['y'] + 1.
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
    result['MAG_AUTO']       = _safe_mag(flux_auto)
    result['MAGERR_AUTO']    = _safe_magerr(flux_auto, fluxerr_auto)
    result['FLUX_PETRO']     = flux_petro
    result['FLUXERR_PETRO']  = fluxerr_petro
    result['MAG_PETRO']      = _safe_mag(flux_petro)
    result['MAGERR_PETRO']   = _safe_magerr(flux_petro, fluxerr_petro)

    # Aperture columns: first aperture has no numeric suffix
    for i, (fl, fle) in enumerate(zip(flux_aper_all, fluxerr_aper_all)):
        suffix = '' if i == 0 else f'_{i}'
        result['FLUX_APER'    + suffix] = fl
        result['FLUXERR_APER' + suffix] = fle
        result['MAG_APER'     + suffix] = _safe_mag(fl)
        result['MAGERR_APER'  + suffix] = _safe_magerr(fl, fle)

    # 12. ASSOC positional matching -------------------------------------------
    if ASSOC_NAME is not None:
        try:
            assoc_data = ascii.read(ASSOC_NAME)
            col_indices = [int(c.strip()) - 1
                           for c in ASSOC_PARAMS.split(',')]
            assoc_x = np.array(assoc_data.columns[col_indices[0]], dtype=float)
            assoc_y = np.array(assoc_data.columns[col_indices[1]], dtype=float)
            # Input positions are 1-indexed; compare with 1-indexed result
            src_x = np.array(result['XWIN_IMAGE'])
            src_y = np.array(result['YWIN_IMAGE'])

            matched_src  = []
            matched_assoc = []
            matched_dist  = []

            for j, (ax, ay) in enumerate(zip(assoc_x, assoc_y)):
                dists   = np.sqrt((src_x - ax)**2 + (src_y - ay)**2)
                nearest = int(np.argmin(dists))
                if dists[nearest] <= ASSOC_RADIUS:
                    matched_src.append(nearest)
                    matched_assoc.append(j)
                    matched_dist.append(float(dists[nearest]))

            if not matched_src:
                if LOGGER:
                    LOGGER.warning(
                        f'ASSOC: no match within {ASSOC_RADIUS:.1f} px in {FITS}')
                return table.Table()

            result = result[matched_src]

            # VECTOR_ASSOC: third column from ASSOC file if available
            if len(col_indices) > 2:
                vec = np.array(assoc_data.columns[col_indices[2]], dtype=float)
                result['VECTOR_ASSOC'] = vec[matched_assoc]
            else:
                result['VECTOR_ASSOC'] = np.zeros(len(matched_src))

            result['NUMBER_ASSOC'] = np.array(matched_assoc, dtype=float) + 1.

        except Exception as exc:
            msg = f'ASSOC matching failed: {exc}'
            print(bcolors.WARNING + msg + bcolors.ENDC)
            if LOGGER:
                LOGGER.warning(msg)
            # Return unfiltered table without ASSOC columns
            result['VECTOR_ASSOC'] = np.zeros(len(result))
            result['NUMBER_ASSOC'] = np.zeros(len(result))

    return result


# ---------------------------------------------------------------------------
# Sextractor output post-processing
# ---------------------------------------------------------------------------

def sextractor_postprocess(DATA, PRINTHELP=True):
    """Post-process source catalog: convert negative fluxes and recompute errors.

    Sets negative fluxes to NaN and recomputes asymmetric magnitude errors
    from flux and flux-error columns.  This makes the catalog consistent with
    forced-photometry upper limits.

    Parameters
    ----------
    DATA : astropy.table.Table
        Source catalog with FLUX_*/FLUXERR_* column pairs.
    PRINTHELP : bool
        Print informational messages about the processing.

    Returns
    -------
    astropy.table.Table
        Modified catalog (in-place modification, also returned).
    """
    if PRINTHELP:
        print(bcolors.WARNING
              + 'Recalculate magnitude errors [ mag_err = -2.5 log(F ± dF)/F ]'
              + bcolors.ENDC)
        print(bcolors.WARNING
              + 'Objects with F ≤ 0 become 3σ upper limits. mag_err = NaN.'
              + bcolors.ENDC)

    flux_keys = [k for k in DATA.colnames
                 if k.startswith('FLUX_') and 'RADIUS' not in k]

    for key in flux_keys:
        err_key = key.replace('FLUX_', 'FLUXERR_')
        if err_key not in DATA.colnames:
            continue
        # Save mask before modifying data (B2 bug fix)
        neg_mask = np.array(DATA[key]) <= 0.
        if np.any(neg_mask):
            DATA[key][neg_mask]     = np.nan
            DATA[err_key][neg_mask] = np.nan

    for key in [k for k in DATA.colnames if k.startswith('FLUXERR_')]:
        flux_key = key.replace('FLUXERR_', 'FLUX_')
        if flux_key not in DATA.colnames:
            continue

        flux    = np.array(DATA[flux_key], dtype=float)
        fluxerr = np.array(DATA[key],      dtype=float)

        with np.errstate(divide='ignore', invalid='ignore'):
            errp = np.where(flux > 0,
                            -2.5 * np.log10(flux - fluxerr)
                            + 2.5 * np.log10(flux),
                            np.nan)
            errm = np.where(flux > 0,
                            -2.5 * np.log10(flux)
                            + 2.5 * np.log10(flux + fluxerr),
                            np.nan)

        mid_key  = key.replace('FLUXERR_', 'MAGERR_')
        errp_key = key.replace('FLUXERR_', 'MAGERRP_')
        errm_key = key.replace('FLUXERR_', 'MAGERRM_')

        DATA[mid_key]  = np.abs(errp + errm) / 2.
        DATA[errp_key] = errp
        DATA[errm_key] = errm

    return DATA


# ---------------------------------------------------------------------------
# Zeropoint determination
# ---------------------------------------------------------------------------

def zeropoint(TABLE_REF, TABLE_NEW, FITS='', LOGGER=None,
              NITER=30000, PATH='', TOLERANCE=1):
    """Compute photometric zeropoints from comparison star measurements.

    Cross-matches the reference catalog (known magnitudes) with the
    instrumental catalog (measured fluxes), then uses bootstrap resampling
    via :func:`stat_tools.statNclip` to derive a robust median zeropoint
    and asymmetric 1σ uncertainties for each aperture.

    Parameters
    ----------
    TABLE_REF : astropy.table.Table
        Reference catalog with columns RA, DEC, MAG_CAT, MAGERR_CAT.
    TABLE_NEW : astropy.table.Table
        Instrumental catalog from :func:`sextractor_photometry` with at least
        ALPHAWIN_J2000 and DELTAWIN_J2000 for cross-matching.
    FITS : str
        Science FITS filename (used for output file naming only).
    LOGGER : logging.Logger, optional
        Logger instance.
    NITER : int
        Number of bootstrap iterations (default: 30 000).
    PATH : str
        Output directory for diagnostic plots.
    TOLERANCE : float
        Cross-matching tolerance in arcseconds (default: 1.0).

    Returns
    -------
    astropy.table.Table
        Zeropoint table with columns METHOD, ZP, ZP_ERRP, ZP_ERRM, NUMBER.
    """
    print(bcolors.OKGREEN
          + f'Bootstrap ZP from {NITER} resamplings\n'
          + bcolors.ENDC)

    # Keep only essential columns in reference table
    for key in list(TABLE_REF.colnames):
        if key not in ('RA', 'DEC', 'MAG_CAT', 'MAGERR_CAT'):
            del TABLE_REF[key]

    TABLE_REF_keys = list(TABLE_REF.colnames)

    # Dynamically build TABLE_NEW column list (do not hardcode aperture count)
    base_keys = [k for k in ('ALPHAWIN_J2000', 'DELTAWIN_J2000',
                              'XWIN_IMAGE', 'YWIN_IMAGE',
                              'MAG_AUTO', 'MAGERR_AUTO',
                              'MAG_PETRO', 'MAGERR_PETRO',
                              'FLUX_AUTO', 'FLUXERR_AUTO',
                              'FLUX_PETRO', 'FLUXERR_PETRO',
                              'FWHM_IMAGE', 'FWHM_WORLD',
                              'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
                              'KRON_RADIUS', 'FLUX_RADIUS', 'FLAGS',
                              'VECTOR_ASSOC', 'NUMBER_ASSOC')
                 if k in TABLE_NEW.colnames]

    aper_keys = sorted(
        [k for k in TABLE_NEW.colnames
         if k.startswith(('MAG_APER', 'MAGERR_APER', 'FLUX_APER', 'FLUXERR_APER'))],
        key=lambda k: (k.split('_')[0],
                       int(k.rsplit('_', 1)[-1]) if k.rsplit('_', 1)[-1].isdigit() else -1))

    TABLE_NEW_keys = base_keys + aper_keys

    # Convert to plain float arrays for fast cross-matching
    TABLE_REF_arr = rfn.structured_to_unstructured(
        np.asarray(TABLE_REF[TABLE_REF_keys]), dtype=float)
    TABLE_NEW_arr = rfn.structured_to_unstructured(
        np.asarray(TABLE_NEW[TABLE_NEW_keys]), dtype=float)

    merged_arr = cat_tools.wrapper_crossmatch(TABLE_REF_arr, TABLE_NEW_arr,
                                               TOLERANCE)
    if len(merged_arr) == 0:
        msg = ('No cross-matched stars for ZP calculation. '
               f'Check coordinates and tolerance ({TOLERANCE} arcsec).')
        print(bcolors.FAIL + msg + bcolors.ENDC)
        if LOGGER:
            LOGGER.error(msg)
        return table.Table(names=('METHOD', 'ZP', 'ZP_ERRP', 'ZP_ERRM', 'NUMBER'),
                           dtype=('U100', 'f', 'f', 'f', 'g'))

    merged = table.Table(merged_arr,
                         names=TABLE_NEW_keys + TABLE_REF_keys + ['DIST'])
    merged['DIST']        = merged['DIST'] * 3600.
    merged['DIST'].format = '.3f'

    keys_mag = [k for k in merged.colnames
                if ('MAG_' in k)
                and k not in ('MAG_CAT', 'MAG_INS')]

    result = table.Table(
        names=('METHOD', 'ZP', 'ZP_ERRP', 'ZP_ERRM', 'NUMBER'),
        dtype=('U100', 'f', 'f', 'f', 'g'))
    result['NUMBER'].format = '7g'

    # Diagnostic ZP plot
    fig = plt.figure(5, figsize=(np.sqrt(2.) * 9, 9))
    fig.subplots_adjust(hspace=0.3, wspace=0.4)

    # Shared axis labels via invisible overlay
    ax_bg = fig.add_subplot(111)
    for spine in ax_bg.spines.values():
        spine.set_visible(False)
    ax_bg.tick_params(labelcolor='w', top=False, bottom=False,
                      left=False, right=False)
    ax_bg.set_xlabel('Apparent magnitude')
    ax_bg.set_ylabel('Zeropoint')
    ax_bg.yaxis.set_label_coords(-0.08, 0.5)

    ncols   = 3
    nrows   = max(1, int(np.ceil(len(keys_mag) / ncols)))

    for i, key in enumerate(keys_mag):
        # Per-star ZP = catalog_mag - instrumental_mag
        temp_zp     = np.array(merged['MAG_CAT'] - merged[key], dtype=float)
        err_key     = 'MAGERR_' + key.split('_', 1)[1]
        if err_key in merged.colnames:
            temp_zp_err = np.sqrt(
                np.array(merged['MAGERR_CAT'], dtype=float)**2
                + np.array(merged[err_key], dtype=float)**2)
        else:
            temp_zp_err = np.array(merged['MAGERR_CAT'], dtype=float)

        pos_mask    = temp_zp > 0
        temp_zp     = np.column_stack([temp_zp[pos_mask],
                                        temp_zp_err[pos_mask]])
        cat_mag_sel = np.array(merged['MAG_CAT'])[pos_mask]

        print(bcolors.OKGREEN + key + bcolors.ENDC)
        print(np.array([f'{v:.3f}' for v in temp_zp[:, 0]]))

        if LOGGER:
            LOGGER.info(key)
            LOGGER.info(np.array([f'{v:.3f}' for v in temp_zp[:, 0]]))

        if len(temp_zp) > 0:
            zp_ana = np.array(stat_tools.statNclip(temp_zp, NITER=NITER))
            result.add_row(np.hstack([key, zp_ana]))
        else:
            result.add_row(np.hstack([key, np.zeros(4)]))

        # Subplot
        ax = fig.add_subplot(nrows, ncols, i + 1)

        if len(temp_zp) > 0:
            zp_25 = np.percentile(temp_zp[:, 0], 25)
            zp_75 = np.percentile(temp_zp[:, 0], 75)
            iqr   = zp_75 - zp_25
            good  = ((temp_zp[:, 0] > zp_25 - 1.5 * iqr)
                     & (temp_zp[:, 0] < zp_75 + 1.5 * iqr))

            ax.axhline(zp_ana[0], lw=4, color=vigit_color_12)
            ax.axhline(zp_ana[0] + zp_ana[1], lw=2, color=vigit_color_12, ls='--')
            ax.axhline(zp_ana[0] - zp_ana[2], lw=2, color=vigit_color_12, ls='--')
            ax.axhline(zp_75 + 1.5 * iqr,     lw=2, color=vigit_color_12, ls=':')
            ax.axhline(zp_25 - 1.5 * iqr,     lw=2, color=vigit_color_12, ls=':')

            if np.any(~good):
                ax.errorbar(cat_mag_sel[~good], temp_zp[~good, 0],
                            temp_zp[~good, 1],
                            marker='o', ms=9, color='0.75',
                            elinewidth=2, capsize=0, lw=0)
            if np.any(good):
                ax.errorbar(cat_mag_sel[good], temp_zp[good, 0],
                            temp_zp[good, 1],
                            marker='o', ms=9, color='k',
                            elinewidth=2, capsize=0, lw=0)

            if len(cat_mag_sel) > 1:
                ax.set_xlim(cat_mag_sel.min() - 0.5, cat_mag_sel.max() + 0.5)
                ax.xaxis.set_major_locator(plt.MultipleLocator(1))

        ax.grid(True)
        label = key.replace('_', r'\_')
        ax.text(0.95, 0.95, label, ha='right', va='top',
                transform=ax.transAxes, fontsize=legend_size)

    plt.savefig(str(Path(PATH) / (Path(FITS).stem + '_zp.pdf')), dpi=600)
    plt.close()

    # FWHM distribution plot
    plt.figure(6, figsize=(np.sqrt(2.) * 9, 9))
    ax_fwhm = plt.subplot(111)

    use_fwhm = ('FWHM_IMAGE' in merged.colnames
                and not np.all(np.array(merged['FWHM_IMAGE']) == 0))
    fwhm_vals = (np.array(merged['FWHM_IMAGE'], dtype=float) if use_fwhm
                 else np.array(merged['FLUX_RADIUS'], dtype=float) * 2. / 1.1)
    fwhm_label = 'FWHM (px)' if use_fwhm else 'FWHM estimate (px)'

    finite = fwhm_vals[np.isfinite(fwhm_vals)]
    if len(finite) > 0:
        ax_fwhm.hist(finite, range=(0, 30), density=True,
                     color=vigit_color_1, bins=20)
        ax_fwhm.hist(finite, range=(0, 30), density=True,
                     color='k', bins=10000, cumulative=True,
                     histtype='step', lw=4)
        p50, p16, p84 = np.percentile(finite, [50, 16, 84])
        ax_fwhm.axvline(p50, lw=4, color='k')
        ax_fwhm.axvline(p16, lw=2, ls=':', color='k')
        ax_fwhm.axvline(p84, lw=2, ls=':', color='k')

    ax_fwhm.set_xlim(0, 30)
    ax_fwhm.set_ylim(0, 1.1)
    ax_fwhm.set_xlabel(fwhm_label)
    ax_fwhm.set_ylabel('Normalised count / CDF')

    plt.savefig(str(Path(PATH) / (Path(FITS).stem + '_fwhm.pdf')), dpi=600)
    plt.close()

    return result
