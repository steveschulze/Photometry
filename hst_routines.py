"""
hst_routines.py — HST-specific aperture photometry routines.

Handles the HST-specific pipeline: curve-of-growth analysis, aperture
corrections via pysynphot, drizzled-image noise corrections, and diagnostic
cutout images.  All other photometry (ground-based) is in extraction.py
and calibration.py.

Public API
----------
hst_aperture_photometry(...)    Aperture photometry on drizzled HST images
hst_make_cutout(...)            Multi-panel aperture diagnostic cutout
hst_cog(...)                    Curve of growth with Monte Carlo errors
hst_scale(instrument, mode)     Native pixel scale per HST detector
hst_zeropoint(header, diameter) AB zeropoint + pysynphot aperture correction
"""

from __future__ import annotations

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"


import sys
from pathlib import Path

import numpy as np
from astropy import stats, table, time
from astropy.io import ascii, fits
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as PathEffects
from scipy import stats as stats_scipy

from extraction import aperture_photometry
from utils import bcolors
from plotsettings import (
    vigit_color_1,
    vigit_color_12,
    color_green,
    color_yellow,
    label_size,
    legend_size,
)


# ---------------------------------------------------------------------------
# HST instrument pixel scales
# ---------------------------------------------------------------------------

def hst_scale(instrument: str, mode: str | None) -> float:
    """Return the native (pre-drizzle) pixel scale for an HST instrument.

    Parameters
    ----------
    instrument : str
        HST instrument name: 'ACS', 'WFC3', 'WFPC2', or 'NICMOS'.
    mode : str | None
        Detector / aperture identifier (e.g. 'WFC', 'UVIS', 'IR').

    Returns
    -------
    float
        Native pixel scale in arcsec/pixel.

    Raises
    ------
    ValueError
        For unknown instrument or mode combinations.
    """
    if instrument == 'NICMOS':
        return 0.05071
    if instrument == 'WFPC2':
        return 0.10

    if instrument == 'ACS':
        if mode and 'WFC' in mode:
            return 0.05
        if mode == 'HRC':
            return 0.0265
        if mode == 'SBC':
            return 0.032
        raise ValueError(f'Unknown ACS mode: {mode!r}')

    if instrument == 'WFC3':
        if mode and 'UVIS' in mode:
            return 0.040
        if mode and 'IR' in mode:
            return 0.13
        raise ValueError(f'Unknown WFC3 mode: {mode!r}')

    raise ValueError(f'Unknown HST instrument: {instrument!r}')


# ---------------------------------------------------------------------------
# HST zeropoint with aperture correction
# ---------------------------------------------------------------------------

def hst_zeropoint(
    header: fits.Header,
    diameter: np.ndarray | list[float],
) -> np.ndarray:
    """Compute HST AB zeropoint with pysynphot aperture correction.

    Derives the infinite-aperture zeropoint from PHOTFLAM and PHOTPLAM,
    then uses pysynphot to compute the enclosed-flux fraction at each
    requested aperture diameter.

    Parameters
    ----------
    header : fits.Header
        FITS header containing HST photometric keywords (PHOTFLAM, PHOTPLAM,
        PHOTMODE or INSTRUME/FILTER/APERTURE/DATE-OBS).
    diameter : array-like
        Aperture diameters in arcseconds.  Diameters ≥ 4″ are treated as
        infinite (no correction).

    Returns
    -------
    np.ndarray
        AB zeropoints (one per diameter) including aperture correction.
    """
    try:
        import pysynphot as pyS
    except ImportError as exc:
        raise ImportError(
            "pysynphot is required for HST zeropoints.\n"
            "Install with: pip install pysynphot"
        ) from exc

    zp_inf = (-2.5 * np.log10(header['PHOTFLAM'])
              - 5.0 * np.log10(header['PHOTPLAM'])
              - 2.408)

    spec_bb  = pyS.BlackBody(10000)
    diameter = np.atleast_1d(
        np.where(np.asarray(diameter) < 4.0, diameter, 4.0))

    try:
        filt = (header['FILTER1']
                if 'CLEAR' in header.get('FILTER2', 'CLEAR')
                else header['FILTER2'])
    except KeyError:
        filt = header.get('FILTER', header.get('FILTNAM1', 'CLEAR'))

    ap_correction = np.ones(len(diameter))

    if header.get('INSTRUME', '').upper() not in ('WFPC2', 'NICMOS'):
        try:
            photmode = header['PHOTMODE'].replace(' ', ',')
            bps     = [pyS.ObsBandpass(f'{photmode},aper#{d:.2f}')
                       for d in diameter]
            bp_ref  = pyS.ObsBandpass(f'{photmode},aper#4.00')
        except Exception:
            mjd  = int(time.Time(header['DATE-OBS'],
                                  format='isot', scale='utc').jd)
            inst = header['INSTRUME']
            det  = header['APERTURE']
            bps  = [pyS.ObsBandpass(
                        f'{inst},{det},{filt},mjd#{mjd},aper#{d:.2f}')
                    for d in diameter]
            bp_ref = pyS.ObsBandpass(
                f'{inst},{det},{filt},mjd#{mjd},aper#4.00')

        cr_ref = pyS.Observation(spec_bb, bp_ref).countrate()
        for i, (d, bp) in enumerate(zip(diameter, bps)):
            if d < 4.0:
                ap_correction[i] = (pyS.Observation(spec_bb, bp).countrate()
                                    / cr_ref)
    else:
        print(bcolors.WARNING
              + f'Aperture correction for {header["INSTRUME"]} not in pysynphot.\n'
              'The AP correction is ≈ −0.1 mag for r = 0.5″. Add by hand.'
              + bcolors.ENDC)

    return zp_inf - (-2.5 * np.log10(ap_correction))


# ---------------------------------------------------------------------------
# HST aperture photometry
# ---------------------------------------------------------------------------

def hst_aperture_photometry(
    fits_path: str | Path,
    positions: np.ndarray,
    radii: np.ndarray,
    inner_annulus: np.ndarray,
    outer_annulus: np.ndarray,
    pix2arcsec: float,
    rms: float,
) -> table.Table:
    """Aperture photometry for HST drizzled images.

    Computes the effective gain from EXPTIME × CCD gain keyword, the
    correlated-noise correction factor for drizzled pixels, and delegates
    to :func:`extraction.aperture_photometry`.

    Parameters
    ----------
    fits_path : str | Path
        HST drizzled FITS file.
    positions : np.ndarray
        (x, y) target position in 1-indexed pixel coordinates.
    radii : np.ndarray
        Aperture radii in arcseconds.
    inner_annulus : np.ndarray
        Inner background annulus radii in arcseconds.
    outer_annulus : np.ndarray
        Outer background annulus radii in arcseconds.
    pix2arcsec : float
        Pixel scale in arcsec/pixel.
    rms : float
        Global image RMS (currently unused, retained for API compatibility).

    Returns
    -------
    astropy.table.Table
        Photometry table (see :func:`extraction.aperture_photometry`).
    """
    with fits.open(str(fits_path)) as hdulist:
        hdu_header = hdulist[0].header.copy()
        if len(hdulist) > 1:
            try:
                hdulist[1].data.shape
                hdu_header += hdulist[1].header
                hdu_data    = hdulist[1].data
            except (AttributeError, TypeError):
                hdu_header += hdulist[1].header
                hdu_data    = hdulist[0].data
        else:
            hdu_data = hdulist[0].data

    radii_px   = np.asarray(radii) / pix2arcsec
    inner_px   = np.asarray(inner_annulus) / pix2arcsec
    outer_px   = np.asarray(outer_annulus) / pix2arcsec

    zp         = hst_zeropoint(hdu_header, 2 * np.asarray(radii))

    for key in ('CCDGAIN', 'ADCGAIN', 'ATODGAIN'):
        if key in hdu_header:
            gain_key = key
            break
    else:
        raise KeyError(f'No gain keyword found in {fits_path}')

    effective_gain = hdu_header['EXPTIME'] * hdu_header[gain_key]

    pixfrac = hdu_header['D001PIXF']
    try:
        native = hst_scale(hdu_header['INSTRUME'], hdu_header['APERTURE'])
    except KeyError:
        native = hst_scale(hdu_header['INSTRUME'], None)

    scale = pix2arcsec / native
    fa    = ((scale / pixfrac) * (1.0 - scale / 3.0 / pixfrac))**2 \
            if scale < pixfrac else (1.0 - pixfrac / 3.0 / scale)**2

    return aperture_photometry(
        hdu_data, positions, radii_px, inner_px, outer_px, rms,
        gain=effective_gain, zeropoint=zp, corr_factor=fa)


# ---------------------------------------------------------------------------
# HST curve of growth
# ---------------------------------------------------------------------------

def hst_cog(
    fits_path: str | Path,
    positions: np.ndarray,
    inner_annulus: float,
    outer_annulus: float,
    pix2arcsec: float,
    rms: float,
    output_dir: str | Path,
) -> None:
    """Compute and plot the curve of growth for an HST source.

    Measures brightness at 20 aperture diameters from 0.2″ to 5″,
    then propagates asymmetric uncertainties via Monte Carlo resampling
    to derive a robust median brightness and scatter.

    Parameters
    ----------
    fits_path : str | Path
        HST drizzled FITS file.
    positions : np.ndarray
        (x, y) target position (1-indexed pixels).
    inner_annulus : float
        Inner annulus scaling factor relative to aperture radius.
    outer_annulus : float
        Outer annulus scaling factor.
    pix2arcsec : float
        Pixel scale in arcsec/pixel.
    rms : float
        Global image RMS.
    output_dir : str | Path
        Output directory for plots and ASCII data.
    """
    apertures_diam = np.linspace(0.2, 5.0, 20)
    inner          = inner_annulus * apertures_diam
    outer          = outer_annulus * apertures_diam

    phot     = hst_aperture_photometry(
        fits_path, positions,
        apertures_diam / 2.0, inner / 2.0, outer / 2.0,
        pix2arcsec, rms)

    mags     = np.array([float(phot[f'MAG_APER_{i}'][0])
                          for i in range(len(apertures_diam))])
    mags_ep  = np.array([float(phot[f'MAGERRP_APER_{i}'][0])
                          for i in range(len(apertures_diam))])
    mags_em  = np.array([float(phot[f'MAGERRM_APER_{i}'][0])
                          for i in range(len(apertures_diam))])

    det      = np.isfinite(mags_ep) & (mags_ep != 0.0)
    ul       = ~det

    out_dir  = Path(output_dir)
    stem     = Path(fits_path).stem

    cog_data = table.Table(
        [apertures_diam, np.round(mags, 3),
         np.round(mags_ep, 3), np.round(mags_em, 3)],
        names=('DIAMETER', 'MAG', 'MAGERRP', 'MAGERRM'))
    ascii.write(cog_data, out_dir / f'{stem}_cog_data.ascii',
                overwrite=True, format='no_header')

    # Monte Carlo error propagation (vectorised)
    niter   = 10_000
    mdet    = mags[det]
    ep_det  = mags_ep[det]
    em_det  = mags_em[det]

    u_mc    = np.random.uniform(size=(niter, len(mdet)))
    values  = np.where(u_mc < 0.5,
                       stats_scipy.norm.ppf(u_mc, mdet, em_det),
                       stats_scipy.norm.ppf(u_mc, mdet, ep_det))
    cog_med = float(np.median(values))
    cog_std = float(np.std(np.median(values, axis=1)))

    cog_stat = table.Table([[cog_med], [cog_med], [cog_std], [niter]],
                            names=('MEAN', 'MEDIAN', 'STD', 'NITER'))
    ascii.write(cog_stat, out_dir / f'{stem}_cog_stat.ascii', overwrite=True)

    # Plot
    fig, ax = plt.subplots(figsize=(9 * np.sqrt(2.), 9))
    for sig, shade in [(3, '0.95'), (2, '0.85'), (1, '0.75')]:
        ax.axhspan(cog_med - sig * cog_std, cog_med + sig * cog_std, color=shade)
    ax.axhline(cog_med, lw=3, color='0.75')

    ax.errorbar(apertures_diam, mags, mags_ep, ms=0, lw=3, color=vigit_color_12)
    ax.errorbar(apertures_diam[det], mags[det],
                [mags_ep[det], mags_em[det]],
                marker='o', ms=9, lw=0, color=vigit_color_12,
                capsize=0, elinewidth=2)
    if ul.any():
        ax.errorbar(apertures_diam[ul], mags[ul],
                    marker='v', ms=9, lw=0, color=vigit_color_12,
                    mec=vigit_color_12)

    ax.text(0.05, 0.95,
            f'median = {cog_med:.3f}, std = {cog_std:.3f}',
            ha='left', va='top', fontsize=label_size,
            transform=ax.transAxes, color='k')
    ax.set_xlabel('Diameter (arcsec)')
    ax.set_ylabel('Brightness (mag, AB)')
    ax.set_xlim(0, apertures_diam[-1])

    valid = mags[np.isfinite(mags)]
    if len(valid):
        ax.set_ylim(valid.max(), valid.min() - 0.5)

    plt.savefig(str(out_dir / f'{stem}_cog.pdf'))
    plt.close()


# ---------------------------------------------------------------------------
# HST diagnostic cutout
# ---------------------------------------------------------------------------

def hst_make_cutout(
    fits_path: str | Path,
    coord_obs: list[float],
    coord_exp: list[float],
    radii: np.ndarray,
    inner_annulus: np.ndarray,
    outer_annulus: np.ndarray,
    pix2arcsec: float,
    output_dir: str | Path,
) -> None:
    """Create a multi-panel aperture diagnostic cutout for HST data.

    Each panel shows the target region at a different aperture size,
    with aperture circle and background annulus overlaid.

    Parameters
    ----------
    fits_path : str | Path
        HST drizzled FITS file.
    coord_obs : list[float]
        Observed (x, y) centroid in 1-indexed pixels.
    coord_exp : list[float]
        Expected (x, y) position from WCS in 1-indexed pixels.
    radii : np.ndarray
        Aperture radii in arcseconds.
    inner_annulus : np.ndarray
        Inner annulus radii in arcseconds.
    outer_annulus : np.ndarray
        Outer annulus radii in arcseconds.
    pix2arcsec : float
        Pixel scale in arcsec/pixel.
    output_dir : str | Path
        Output directory for the PDF plot.
    """
    with fits.open(str(fits_path)) as hdulist:
        if len(hdulist) > 1:
            try:
                hdulist[1].data.shape
                image = hdulist[1].data
            except (AttributeError, TypeError):
                image = hdulist[0].data
        else:
            image = hdulist[0].data

    halfwidth = 100
    cx, cy    = int(coord_obs[0]), int(coord_obs[1])
    h, w      = image.shape
    xmin, xmax = max(cx - halfwidth, 0), min(cx + halfwidth, w)
    ymin, ymax = max(cy - halfwidth, 0), min(cy + halfwidth, h)
    cutout     = image[ymin:ymax, xmin:xmax]

    flat = cutout.flatten()
    flat = flat[np.isfinite(flat) & (flat > 0)]
    vmin = np.percentile(flat, 55) if len(flat) else 1e-3
    vmax = np.percentile(flat, 99) if len(flat) else 1.0

    n    = len(radii)
    fig, axes = plt.subplots(
        max(1, int(np.ceil(n / 3))), 3,
        figsize=(np.sqrt(2.) * 9, 9))
    axes_flat = np.array(axes).flatten()

    for i in range(n):
        ax = axes_flat[i]
        ax.imshow(cutout, cmap='binary', interpolation='nearest',
                  origin='lower', norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.plot(halfwidth, halfwidth, marker='x', mew=5,
                color=color_green, ms=12)

        if len(coord_obs) > 0:
            ox = coord_obs[0] - cx + halfwidth - 1
            oy = coord_obs[1] - cy + halfwidth - 1
            ax.errorbar(ox, oy, mew=5, marker='x', color=vigit_color_1, ms=12)

        for r, col in [
            (radii[i] / pix2arcsec,         vigit_color_1),
            (inner_annulus[i] / pix2arcsec,  color_yellow),
            (outer_annulus[i] / pix2arcsec,  color_yellow),
        ]:
            ax.add_patch(plt.Circle((halfwidth, halfwidth), r,
                                     ls='solid', lw=5, fc='None', ec=col))

        ax.set_xlim(0, 2 * halfwidth)
        ax.set_ylim(0, 2 * halfwidth)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)

        ax.text(0.05, 0.95,
                f"\\textbf{{Diam:}} $\\mathbf{{ {2*radii[i]:.1f}''}}$",
                ha='left', va='top', transform=ax.transAxes,
                fontsize=legend_size - 4,
                path_effects=[PathEffects.withStroke(linewidth=6, foreground='w')])

    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.savefig(
        str(Path(output_dir) / (Path(fits_path).stem + '.pdf')), dpi=600)
    plt.close()
