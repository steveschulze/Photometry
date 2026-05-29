"""
calibration.py — Photometric calibration, zeropoint determination, and output.

Builds the local comparison-star sequence, determines the photometric
zeropoint via bootstrap resampling, and assembles the final science and
diagnostic output tables and plots.

Public API
----------
local_sequence(catalog, ...)     Select and write comparison stars
zeropoint(ref_table, phot_table, ...) Bootstrap ZP per aperture
make_scicat(...)                 Assemble the science photometry output table
make_poststamp(...)              Two-panel diagnostic cutout image
"""

from __future__ import annotations

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"


import logging
import sys
from pathlib import Path

import numpy as np
import numpy.lib.recfunctions as rfn
from astropy import stats, table, time
from astropy.io import ascii, fits
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as PathEffects
from scipy import optimize

import cat_tools
import fits_tools
from utils import bcolors, bootstrap_zp_stats
from plotsettings import (
    vigit_color_1,
    vigit_color_12,
    color_green,
    label_size,
    legend_size,
    left, right, top, bottom,
)


# ---------------------------------------------------------------------------
# Local sequence selection
# ---------------------------------------------------------------------------

def local_sequence(
    catalog: table.Table,
    *,
    auto: bool = False,
    output_file: str | Path | None = None,
    fits_path: str | Path = '',
    logger: logging.Logger | None = None,
    lower: float = 10,
    output_dir: str | Path = '',
    upper: float = 90,
) -> dict[str, int | table.Table]:
    """Select comparison stars for the local photometric sequence.

    In *auto* mode the magnitude range is read from an instrument/filter
    lookup table derived from the FITS filename.  In interactive mode the
    user is prompted to apply custom brightness cuts.

    Parameters
    ----------
    catalog : astropy.table.Table
        Cross-matched catalog.  Must contain columns ``MAG_INS``, ``MAG_CAT``,
        ``MAGERR_CAT``, ``MAGERR_INS``, ``XWIN_IMAGE``, ``YWIN_IMAGE``,
        ``ALPHAWIN_J2000``, ``DELTAWIN_J2000``.
    auto : bool
        Use the automatic magnitude range from the instrument lookup table.
    output_file : str | Path | None
        Path for the cleaned local-sequence ASCII catalog.
    fits_path : str | Path
        Science FITS filename — parsed to identify telescope/filter in auto mode.
    logger : logging.Logger | None
        Logger instance.
    lower : float
        Lower percentile cut on instrumental magnitude (manual mode).
    output_dir : str | Path
        Directory for diagnostic plots.
    upper : float
        Upper percentile cut on instrumental magnitude (manual mode).

    Returns
    -------
    dict with keys:
        ``'NUMSTARS'``: int — number of stars in the cleaned sequence.
        ``'CAT'``: astropy.table.Table — cleaned comparison-star catalog.
    """
    if output_file is None:
        print(bcolors.FAIL + 'output_file not specified for local_sequence.'
              + bcolors.ENDC)
        sys.exit(1)

    output_file = Path(output_file)
    fits_path   = Path(fits_path)
    out_dir     = Path(output_dir)

    # Instrument/filter magnitude limits (used in auto mode)
    _limits = {
        ('PanSTARRS', 'g'): (-18, -13), ('PanSTARRS', 'r'): (-18, -13),
        ('PanSTARRS', 'i'): (-18, -13), ('PanSTARRS', 'z'): (-18, -12),
        ('PanSTARRS', 'y'): (-18, -12),
        ('2MASS', 'J'): (-10, -6),  ('2MASS', 'H'): (-10, -6),
        ('2MASS', 'K'): (-10, -6.5),
        ('UKIDSS', 'J'): (-16, -11), ('UKIDSS', 'H'): (-16, -11),
        ('UKIDSS', 'K'): (-16, -11), ('UKIDSS', 'Y'): (-16, -8),
        ('SDSS', 'u'): (-7, -3),  ('SDSS', 'g'): (-10, -3),
        ('SDSS', 'r'): (-10, -3), ('SDSS', 'i'): (-10, -4),
        ('SDSS', 'z'): (-10, -4),
    }

    def _fit_offset(x: np.ndarray, y: np.ndarray,
                    yerr: np.ndarray) -> float:
        """Fit a constant offset  y = offset + x."""
        fitfunc = lambda p, xi: p + xi
        if np.all(yerr != 0):
            out = optimize.leastsq(
                lambda p, xi, yi, e: (yi - fitfunc(p, xi)) / e,
                0.0, args=(x, y, yerr), full_output=True)
        else:
            out = optimize.leastsq(
                lambda p, xi, yi: yi - fitfunc(p, xi),
                0.0, args=(x, y), full_output=True)
        return float(out[0][0])

    def _save_cat(cat_out: table.Table) -> None:
        """Write the local-sequence catalog to *output_file* in ASCII format."""
        ascii.write(
            np.column_stack([
                cat_out['XWIN_IMAGE'], cat_out['YWIN_IMAGE'],
                cat_out['ALPHAWIN_J2000'], cat_out['DELTAWIN_J2000'],
                cat_out['MAG_CAT'], cat_out['MAGERR_CAT'],
            ]),
            str(output_file),
            names=['XWIN_IMAGE', 'YWIN_IMAGE',
                   'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'MAG', 'MAG_ERR'],
            overwrite=True)

    def _diagnostic_plot(cat_all: table.Table, cat_sel: table.Table,
                         zp_off: float,
                         bright_cut: float | None, faint_cut: float | None,
                         fignum: int, show: bool = False) -> None:
        """Plot instrumental vs. apparent magnitude for the local sequence.

        Highlights selected stars (black) against all detected candidates
        (grey), with the best-fit constant offset and optional bright/faint
        magnitude cut lines overlaid.
        """
        plt.figure(fignum, figsize=(9 * np.sqrt(2.), 9))
        ax = plt.subplot(111)

        mag_ins = np.asarray(cat_all['MAG_INS'])
        xlim    = np.array([mag_ins.min() - 0.2, mag_ins.min() + 0.2])
        xlim[1] = max(xlim[1], mag_ins.max() + 0.2)

        ax.fill_between(xlim, xlim + zp_off - 0.1, xlim + zp_off + 0.1,
                        color=vigit_color_12, alpha=0.25)
        ax.plot(xlim, xlim + zp_off, lw=2, color=vigit_color_12)

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

        ax.set_xlabel("Instrumental magnitude")
        ax.set_ylabel("Apparent magnitude")
        ax.grid(True)

        if len(cat_sel['MAG_CAT']) > 1:
            ax.set_xlim(mag_ins.min() - 0.2,
                        min(mag_ins.max(), 0) + 0.2)
            ax.set_ylim(np.asarray(cat_sel['MAG_CAT']).min() - 0.2,
                        np.asarray(cat_sel['MAG_CAT']).max() + 0.2)

        plt.savefig(str(out_dir / (fits_path.stem + '_std.pdf')), dpi=600)
        if show:
            plt.show()
        plt.close()

    # ---- AUTO mode -----------------------------------------------------------
    if auto:
        parts     = fits_path.stem.split('_')
        tel       = parts[1] if len(parts) > 1 else ''
        filt_list = [p for p in parts[2:] if len(p) == 1]
        filt      = filt_list[0] if filt_list else ''

        limits    = _limits.get((tel, filt))
        if limits is None:
            raise ValueError(
                f'No magnitude limits defined for telescope={tel!r}, '
                f'filter={filt!r}')
        bright, faint = limits

        print(bcolors.OKGREEN
              + f'AUTO magnitude range: {bright:.2f} to {faint:.2f}'
              + bcolors.ENDC)
        if logger:
            logger.info(f'Magnitude range: {bright:.2f} - {faint:.2f}')

        mask    = ((catalog['MAG_INS'] >= bright)
                   & (catalog['MAG_INS'] <= faint))
        cat_sel = catalog[mask]

        zp_off = _fit_offset(np.asarray(cat_sel['MAG_INS']),
                              np.asarray(cat_sel['MAG_CAT']),
                              np.asarray(cat_sel['MAGERR_CAT']))
        _diagnostic_plot(catalog, cat_sel, zp_off, bright, faint, fignum=2)

        print(bcolors.OKGREEN
              + f'Number of stars: {mask.sum()}'
              + bcolors.ENDC)
        _save_cat(cat_sel)
        return {'NUMSTARS': int(mask.sum()), 'CAT': cat_sel}

    # ---- Full range (no interactive cut) ------------------------------------
    if lower == 0.0 and upper == 100:
        zp_off = _fit_offset(np.asarray(catalog['MAG_INS']),
                              np.asarray(catalog['MAG_CAT']),
                              np.asarray(catalog['MAGERR_CAT']))
        _diagnostic_plot(catalog, catalog, zp_off, None, None, fignum=2)
        _save_cat(catalog)
        return {'NUMSTARS': len(catalog), 'CAT': catalog}

    # ---- Interactive percentile cut -----------------------------------------
    mag_ins = np.asarray(catalog['MAG_INS'])
    if len(mag_ins) > 10:
        bright_cut = float(np.percentile(mag_ins, lower))
        faint_cut  = float(np.percentile(mag_ins, upper))
    else:
        bright_cut = float(mag_ins.min())
        faint_cut  = float(mag_ins.max())

    mask    = (mag_ins >= bright_cut) & (mag_ins <= faint_cut)
    cat_sel = catalog[mask]
    zp_off  = _fit_offset(np.asarray(cat_sel['MAG_INS']),
                           np.asarray(cat_sel['MAG_CAT']),
                           np.asarray(cat_sel['MAGERR_CAT']))

    _diagnostic_plot(catalog, cat_sel, zp_off, bright_cut, faint_cut,
                     fignum=2, show=True)

    print(f'Current cuts:  bright = {bright_cut:.2f},  faint = {faint_cut:.2f}')
    var1 = input(bcolors.BOLD + bcolors.OKBLUE
                 + '\nApply a magnitude cut? [y|[n]] ' + bcolors.ENDC)

    cat_cleaned = table.Table()

    if var1 != 'y':
        _save_cat(cat_sel)
    else:
        while var1 == 'y':
            bright_cut = float(input(bcolors.WARNING + 'Bright limit: ' + bcolors.ENDC))
            faint_cut  = float(input(bcolors.WARNING + 'Faint limit:  ' + bcolors.ENDC))

            cat_cleaned = catalog[
                (mag_ins >= bright_cut) & (mag_ins <= faint_cut)]
            print(bcolors.OKGREEN
                  + f'\nStars: {len(cat_cleaned)}\n' + bcolors.ENDC)

            if len(cat_cleaned):
                zp_off = _fit_offset(
                    np.asarray(cat_cleaned['MAG_INS']),
                    np.asarray(cat_cleaned['MAG_CAT']),
                    np.asarray(cat_cleaned['MAGERR_CAT']))
            _diagnostic_plot(catalog, cat_cleaned, zp_off,
                              bright_cut, faint_cut, fignum=3, show=True)
            _save_cat(cat_cleaned)

            var1 = input(bcolors.BOLD + bcolors.WARNING
                         + 'Change cuts? [y|[n]] ' + bcolors.ENDC)

    if logger:
        logger.info(f'Magnitude range: {bright_cut:.2f} - {faint_cut:.2f}')

    final = cat_cleaned if len(cat_cleaned) else cat_sel
    if len(final) == 0:
        print(bcolors.FAIL + 'Empty catalog after cuts.' + bcolors.ENDC)
        sys.exit(1)

    return {'NUMSTARS': len(final), 'CAT': final}


# ---------------------------------------------------------------------------
# Zeropoint determination
# ---------------------------------------------------------------------------

def zeropoint(
    ref_table: table.Table,
    phot_table: table.Table,
    *,
    fits_path: str | Path = '',
    logger: logging.Logger | None = None,
    niter: int = 30_000,
    output_dir: str | Path = '',
    tolerance: float = 1.0,
) -> table.Table:
    """Determine photometric zeropoints from comparison star measurements.

    Cross-matches the reference catalog (known magnitudes) against the
    instrumental catalog (measured fluxes), then runs
    :func:`utils.bootstrap_zp_stats` for each aperture.

    Parameters
    ----------
    ref_table : astropy.table.Table
        Reference catalog — columns ``RA``, ``DEC``, ``MAG_CAT``, ``MAGERR_CAT``.
    phot_table : astropy.table.Table
        Instrumental photometry from :func:`extraction.extract_sources` with at
        least ``ALPHAWIN_J2000`` and ``DELTAWIN_J2000``.
    fits_path : str | Path
        Science FITS filename — used only for output file naming.
    logger : logging.Logger | None
        Logger instance.
    niter : int
        Bootstrap iterations (default: 30 000).
    output_dir : str | Path
        Directory for diagnostic PDF plots.
    tolerance : float
        Cross-matching tolerance in arcseconds (default: 1.0).

    Returns
    -------
    astropy.table.Table
        Columns: ``METHOD``, ``ZP``, ``ZP_ERRP``, ``ZP_ERRM``, ``NUMBER``.
    """
    fits_path = Path(fits_path)
    out_dir   = Path(output_dir)

    print(bcolors.OKGREEN + f'Bootstrap ZP from {niter} resamplings\n' + bcolors.ENDC)

    # Keep only essential reference columns
    for key in list(ref_table.colnames):
        if key not in ('RA', 'DEC', 'MAG_CAT', 'MAGERR_CAT'):
            del ref_table[key]

    ref_keys = list(ref_table.colnames)

    # Build phot_table column list dynamically (no hardcoded aperture count)
    _base = ['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE',
             'MAG_AUTO', 'MAGERR_AUTO', 'MAG_PETRO', 'MAGERR_PETRO',
             'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_PETRO', 'FLUXERR_PETRO',
             'FWHM_IMAGE', 'FWHM_WORLD', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
             'KRON_RADIUS', 'FLUX_RADIUS', 'FLAGS', 'VECTOR_ASSOC',
             'NUMBER_ASSOC']
    base_keys = [k for k in _base if k in phot_table.colnames]

    aper_keys = sorted(
        [k for k in phot_table.colnames
         if k.startswith(('MAG_APER', 'MAGERR_APER', 'FLUX_APER', 'FLUXERR_APER'))],
        key=lambda k: (
            k.split('_')[0],
            int(k.rsplit('_', 1)[-1]) if k.rsplit('_', 1)[-1].isdigit() else -1
        ))

    phot_keys = base_keys + aper_keys

    # Convert to float arrays for crossmatch
    ref_arr  = rfn.structured_to_unstructured(
        np.asarray(ref_table[ref_keys]), dtype=float)
    phot_arr = rfn.structured_to_unstructured(
        np.asarray(phot_table[phot_keys]), dtype=float)

    merged_arr = cat_tools.crossmatch_catalogs(ref_arr, phot_arr, tolerance)

    if len(merged_arr) == 0:
        msg = (f'No cross-matches for ZP calculation '
               f'(tolerance = {tolerance} arcsec). '
               'Check image WCS and reference catalogue.')
        print(bcolors.FAIL + msg + bcolors.ENDC)
        if logger:
            logger.error(msg)
        return table.Table(names=('METHOD', 'ZP', 'ZP_ERRP', 'ZP_ERRM', 'NUMBER'),
                           dtype=('U100', 'f', 'f', 'f', 'g'))

    merged          = table.Table(merged_arr, names=phot_keys + ref_keys + ['DIST'])
    merged['DIST']  = merged['DIST'] * 3600.0
    merged['DIST'].format = '.3f'

    mag_keys = [k for k in merged.colnames
                if ('MAG_' in k) and k not in ('MAG_CAT', 'MAG_INS')]

    result = table.Table(
        names=('METHOD', 'ZP', 'ZP_ERRP', 'ZP_ERRM', 'NUMBER'),
        dtype=('U100', 'f', 'f', 'f', 'g'))
    result['NUMBER'].format = '7g'

    # ZP diagnostic figure
    fig    = plt.figure(5, figsize=(np.sqrt(2.) * 9, 9))
    fig.subplots_adjust(hspace=0.3, wspace=0.4)
    ax_bg  = fig.add_subplot(111)
    for spine in ax_bg.spines.values():
        spine.set_visible(False)
    ax_bg.tick_params(labelcolor='w', top=False, bottom=False,
                      left=False, right=False)
    ax_bg.set_xlabel('Apparent magnitude (mag)')
    ax_bg.set_ylabel(r'$\mathrm{ZP} = m_\mathrm{cat} - m_\mathrm{inst}$ (mag)')
    ax_bg.yaxis.set_label_coords(-0.08, 0.5)

    ncols = 3
    nrows = max(1, int(np.ceil(len(mag_keys) / ncols)))

    for i, key in enumerate(mag_keys):
        # ZP = m_cat - m_inst  (positive; e.g. ~22.3 for a typical ground frame)
        zp_vals = np.asarray(merged['MAG_CAT'] - merged[key], dtype=float)
        err_key = 'MAGERR_' + key.split('_', 1)[1]
        if err_key in merged.colnames:
            zp_errs = np.sqrt(
                np.asarray(merged['MAGERR_CAT'], dtype=float)**2
                + np.asarray(merged[err_key], dtype=float)**2)
        else:
            zp_errs = np.asarray(merged['MAGERR_CAT'], dtype=float)

        pos     = zp_vals > 0
        zp_arr  = np.column_stack([zp_vals[pos], zp_errs[pos]])
        cat_sel = np.asarray(merged['MAG_CAT'])[pos]

        print(bcolors.OKGREEN + key + bcolors.ENDC)
        print(np.array([f'{v:.3f}' for v in zp_arr[:, 0]]))
        if logger:
            logger.info('%s: %s', key, [f'{v:.3f}' for v in zp_arr[:, 0]])

        if len(zp_arr) > 0:
            zp_ana = bootstrap_zp_stats(zp_arr, niter=niter)
            result.add_row(np.hstack([key, zp_ana]))
        else:
            result.add_row(np.hstack([key, np.zeros(4)]))

        ax = fig.add_subplot(nrows, ncols, i + 1)

        if len(zp_arr) > 0:
            p25, p75 = np.percentile(zp_arr[:, 0], [25, 75])
            iqr      = p75 - p25
            good     = ((zp_arr[:, 0] > p25 - 1.5 * iqr)
                        & (zp_arr[:, 0] < p75 + 1.5 * iqr))

            ax.axhline(zp_ana[0],              lw=4, color=vigit_color_12)
            ax.axhline(zp_ana[0] + zp_ana[1], lw=2, color=vigit_color_12, ls='--')
            ax.axhline(zp_ana[0] - zp_ana[2], lw=2, color=vigit_color_12, ls='--')
            ax.axhline(p75 + 1.5 * iqr,       lw=2, color=vigit_color_12, ls=':')
            ax.axhline(p25 - 1.5 * iqr,       lw=2, color=vigit_color_12, ls=':')

            if np.any(~good):
                ax.errorbar(cat_sel[~good], zp_arr[~good, 0], zp_arr[~good, 1],
                            marker='o', ms=9, color='0.75',
                            elinewidth=2, capsize=0, lw=0)
            if np.any(good):
                ax.errorbar(cat_sel[good], zp_arr[good, 0], zp_arr[good, 1],
                            marker='o', ms=9, color='k',
                            elinewidth=2, capsize=0, lw=0)
            if len(cat_sel) > 1:
                ax.set_xlim(cat_sel.min() - 0.5, cat_sel.max() + 0.5)
                ax.xaxis.set_major_locator(plt.MultipleLocator(1))

        ax.grid(True)
        ax.text(0.95, 0.95, key.replace('_', r'\_'),
                ha='right', va='top', transform=ax.transAxes,
                fontsize=legend_size)

    plt.savefig(str(out_dir / (fits_path.stem + '_zp.pdf')), dpi=600)
    plt.close()

    # FWHM distribution plot
    plt.figure(6, figsize=(np.sqrt(2.) * 9, 9))
    ax_fwhm = plt.subplot(111)

    use_fwhm  = ('FWHM_IMAGE' in merged.colnames
                 and not np.all(np.asarray(merged['FWHM_IMAGE']) == 0))
    fwhm_vals = (np.asarray(merged['FWHM_IMAGE'], dtype=float) if use_fwhm
                 else np.asarray(merged['FLUX_RADIUS'], dtype=float) * 2.0 / 1.1)
    finite = fwhm_vals[np.isfinite(fwhm_vals)]

    if len(finite):
        ax_fwhm.hist(finite, range=(0, 30), density=True, bins=20,
                     color=vigit_color_1)
        ax_fwhm.hist(finite, range=(0, 30), density=True, bins=10000,
                     cumulative=True, histtype='step', lw=4, color='k')
        p16, p50, p84 = np.percentile(finite, [16, 50, 84])
        ax_fwhm.axvline(p50, lw=4, color='k')
        ax_fwhm.axvline(p16, lw=2, ls=':', color='k')
        ax_fwhm.axvline(p84, lw=2, ls=':', color='k')

    ax_fwhm.set_xlim(0, 30)
    ax_fwhm.set_ylim(0, 1.1)
    ax_fwhm.set_xlabel('FWHM (px)' if use_fwhm else 'FWHM estimate (px)')
    ax_fwhm.set_ylabel('Normalised count / CDF')

    plt.savefig(str(out_dir / (fits_path.stem + '_fwhm.pdf')), dpi=600)
    plt.close()

    return result


# ---------------------------------------------------------------------------
# Science output catalog builder
# ---------------------------------------------------------------------------

def make_scicat(
    fits_path: str | Path,
    object_props: table.Table,
    source_phot: table.Table,
    forced_phot: table.Table | list,
    zp_summary: table.Table,
    search_radius: float,
    logger: logging.Logger | None = None,
) -> table.Table:
    """Assemble the final science photometry summary table.

    Parameters
    ----------
    fits_path : str | Path
        Science FITS filename.
    object_props : astropy.table.Table
        Target properties: RA, DEC, X_EXP, Y_EXP, OBJECT, PHOTCAL.
    source_phot : astropy.table.Table
        Detections near the target (empty table if non-detected).
    forced_phot : astropy.table.Table or list
        Forced photometry at the target position (used when source_phot is empty).
    zp_summary : astropy.table.Table
        Zeropoint table from :func:`zeropoint`.
    search_radius : float
        Detection search radius in arcseconds.
    logger : logging.Logger | None
        Logger instance.

    Returns
    -------
    astropy.table.Table
        Summary table with columns PROPERTY, VALUE, ERROR+, ERROR−, COMMENT.
    """
    hdu_header = fits.getheader(str(fits_path))
    catalog    = table.Table(
        names=('PROPERTY', 'VALUE', 'ERROR+', 'ERROR-', 'COMMENT'),
        dtype=('U100', 'f', 'f', 'f', 'U100'))

    catalog.add_row(['FILENAME', np.nan, np.nan, np.nan, str(fits_path)])

    for key in ('DATE-OBS', 'MJD', 'EXPTIME', 'NCOMBINE'):
        if key in hdu_header:
            catalog.add_row([key, np.nan, np.nan, np.nan, str(hdu_header[key])])
        else:
            catalog.add_row([key, np.nan, np.nan, np.nan,
                             '1' if key == 'NCOMBINE' else '...'])

    if 'OBSMJD' in hdu_header:
        catalog['COMMENT'][catalog['PROPERTY'] == 'DATE-OBS'] = \
            time.Time(hdu_header['OBSMJD'], format='mjd').isot
        catalog['COMMENT'][catalog['PROPERTY'] == 'MJD'] = str(hdu_header['OBSMJD'])

    # Attempt to recover DATE-OBS from non-standard keywords
    if 'T' not in str(catalog['COMMENT'][catalog['PROPERTY'] == 'DATE-OBS'][0]):
        for dk, tk in (('ut-date', 'ut-time'), ('UT-DATE', 'UT-TIME')):
            if dk in hdu_header and tk in hdu_header:
                try:
                    utc = f"{hdu_header[dk]}T{hdu_header[tk]}"
                    catalog['COMMENT'][catalog['PROPERTY'] == 'DATE-OBS'] = utc
                    mjd = time.Time(utc, format='isot', scale='utc').mjd
                    catalog['COMMENT'][catalog['PROPERTY'] == 'MJD'] = f'{mjd:.7f}'
                    break
                except Exception:
                    pass

    catalog.add_row(['PHOTCAL', np.nan, np.nan, np.nan,
                     str(object_props['PHOTCAL'][0])])
    catalog.add_row(['RA',  float(f'{object_props["RA"][0]:.6f}'),
                     np.nan, np.nan, 'degree'])
    catalog.add_row(['DEC', float(f'{object_props["DEC"][0]:.6f}'),
                     np.nan, np.nan, 'degree'])
    catalog.add_row(['XWIN_IMAGE_EXP', float(f'{object_props["X_EXP"][0]:.1f}'),
                     np.nan, np.nan, 'pixel'])
    catalog.add_row(['YWIN_IMAGE_EXP', float(f'{object_props["Y_EXP"][0]:.1f}'),
                     np.nan, np.nan, 'pixel'])

    pix_scale = fits_tools.pix2arcsec(str(fits_path))

    if len(source_phot) > 0:
        msg = (f'One or more objects found within {search_radius:.1f} arcsec '
               f'from the source position')
        print(bcolors.OKGREEN + f'\n{msg}\n' + bcolors.ENDC)
        if logger:
            logger.info(msg)

        x_exp = float(object_props['X_EXP'][0])
        y_exp = float(object_props['Y_EXP'][0])
        dist  = [float(np.hypot(float(r['XWIN_IMAGE']) - x_exp,
                                 float(r['YWIN_IMAGE']) - y_exp))
                 for r in source_phot]

        source_phot['DISTANCE (px)']     = dist
        source_phot['DISTANCE (arcsec)'] = (
            source_phot['DISTANCE (px)'] * pix_scale)
        source_phot.sort('DISTANCE (px)')

        for key in ('XWIN_IMAGE', 'YWIN_IMAGE'):
            catalog.add_row([f'{key}_OBS',
                             float(f'{source_phot[key][0]:.1f}'),
                             np.nan, np.nan, 'pixel'])
        for key in ('ALPHAWIN_J2000', 'DELTAWIN_J2000'):
            catalog.add_row([f'{key}_OBS',
                             float(f'{source_phot[key][0]:.5f}'),
                             np.nan, np.nan, 'degree'])

        catalog.add_row(['DISTANCE (px)',
                         float(f'{source_phot["DISTANCE (px)"][0]:.2f}'),
                         np.nan, np.nan, 'pixel'])
        catalog.add_row(['DISTANCE (arcsec)',
                         float(f'{source_phot["DISTANCE (arcsec)"][0]:.2f}'),
                         np.nan, np.nan, 'arcsec'])
        catalog.add_row(['FLUX_RADIUS (arcsec)',
                         float(f'{source_phot["FLUX_RADIUS"][0] * pix_scale:.2f}'),
                         np.nan, np.nan, 'arcsec'])

        for key in [k for k in source_phot.colnames if 'MAG_' in k]:
            errp_key = key.replace('MAG_', 'MAGERRP_')
            errm_key = key.replace('MAG_', 'MAGERRM_')
            catalog.add_row([
                key,
                float(f'{source_phot[key][0]:.3f}'),
                float(f'{source_phot[errp_key][0]:.3f}'),
                float(f'{source_phot[errm_key][0]:.3f}'),
                'mag'])

    else:
        msg = (f'No object found within {search_radius:.1f} arcsec '
               f'from the source position')
        print(bcolors.WARNING + f'\n{msg}\n' + bcolors.ENDC)
        if logger:
            logger.info(msg)

        for key in ('XWIN_IMAGE', 'YWIN_IMAGE'):
            catalog.add_row([f'{key}_OBS', np.nan, np.nan, np.nan, 'pixel'])
        catalog.add_row(['DISTANCE (px)',     np.nan, np.nan, np.nan, 'pixel'])
        catalog.add_row(['DISTANCE (arcsec)', np.nan, np.nan, np.nan, 'arcsec'])

        # Forced photometry — iterate over apertures
        if not isinstance(forced_phot, list) and len(forced_phot) > 0:
            mag_keys = sorted(
                [k for k in forced_phot.colnames if k.startswith('MAG_APER_')],
                key=lambda k: int(k.split('_')[-1]))
            for i, key in enumerate(mag_keys):
                idx      = key.split('_')[-1]
                fnuerr   = float(forced_phot[f'FNUERR_APER_{idx}'][0])
                mag_3sig = 23.9 - 2.5 * np.log10(3 * fnuerr)
                suffix   = '' if i == 0 else str(i)

                catalog.add_row([f'MAG_APER{suffix}_PHOTUTILS',
                                  float(f'{forced_phot[key][0]:.3f}'),
                                  float(f'{forced_phot[f"MAGERRP_APER_{idx}"][0]:.3f}'),
                                  float(f'{forced_phot[f"MAGERRM_APER_{idx}"][0]:.3f}'),
                                  'mag'])
                catalog.add_row([f'MAG_APER{suffix}_PHOTUTILS_3SIGMA',
                                  float(f'{mag_3sig:.3f}'),
                                  np.nan, np.nan, 'mag'])
                catalog.add_row([f'FNU_APER{suffix}_PHOTUTILS',
                                  float(f'{forced_phot[f"FNU_APER_{idx}"][0]:.3e}'),
                                  fnuerr, fnuerr, 'microJy'])

    return catalog


# ---------------------------------------------------------------------------
# Diagnostic poststamp
# ---------------------------------------------------------------------------

def make_poststamp(
    fits_path: str | Path,
    coord_exp: tuple[float, float],
    coord_obs: tuple[float | None, float | None],
    output_dir: str | Path = '',
) -> None:
    """Create a two-panel diagnostic poststamp.

    Shows the background-subtracted check image (left panel) and the science
    image (right panel) centred on the expected target position, with crosses
    marking the expected (blue) and observed (green) positions.

    Parameters
    ----------
    fits_path : str | Path
        Science FITS filename.
    coord_exp : tuple[float, float]
        Expected (x, y) position from WCS (1-indexed pixels).
    coord_obs : tuple[float | None, float | None]
        Observed centroid from source detection (1-indexed pixels).
        Pass ``(None, None)`` if the source was not detected.
    output_dir : str | Path
        Directory for the output PDF (default: same as FITS file).
    """
    fits_path  = Path(fits_path)
    out_dir    = Path(output_dir) if output_dir else fits_path.parent
    check_path = out_dir / ('check_' + fits_path.name)
    halfwidth  = 50
    cx, cy     = int(coord_exp[0]), int(coord_exp[1])

    def _load(path: Path) -> tuple[np.ndarray, int, int]:
        with fits.open(str(path)) as h:
            img = h[1].data if (len(h) > 1 and h[1].data is not None) else h[0].data
            hdr = h[1].header if len(h) > 1 else h[0].header
        return img, int(hdr['NAXIS1']), int(hdr['NAXIS2'])

    def _vrange(img: np.ndarray, cx: int, cy: int) -> tuple[float, float]:
        ny, nx = img.shape
        cut    = img[max(cy - halfwidth, 0):min(cy + halfwidth, ny),
                     max(cx - halfwidth, 0):min(cx + halfwidth, nx)]
        flat   = np.asarray(sorted(cut.flatten()[np.isfinite(cut.flatten())]))
        if not len(flat):
            return 1e-3, 1.0
        pcts = np.array([np.percentile(flat, p) for p in range(50, 99)])
        vmin = pcts[pcts > 0][1] if np.any(pcts > 0) else float(flat.max()) * 2
        return float(vmin), float(np.percentile(flat, 99))

    check_img, nx, ny = _load(check_path if check_path.exists() else fits_path)
    sci_img, _, _     = _load(fits_path)

    def _cut(img: np.ndarray) -> np.ndarray:
        h, w = img.shape
        return img[max(cy - halfwidth, 0):min(cy + halfwidth, h),
                   max(cx - halfwidth, 0):min(cx + halfwidth, w)]

    check_cut = _cut(check_img)
    sci_cut   = _cut(sci_img)

    fig, (ax_chk, ax_sci) = plt.subplots(1, 2, figsize=(np.sqrt(2.) * 9, 9))

    for ax, img_cut, src_img in [
        (ax_chk, check_cut, check_img),
        (ax_sci, sci_cut,   sci_img),
    ]:
        vmin, vmax = _vrange(src_img, cx, cy)
        ax.imshow(img_cut, cmap='binary', interpolation='nearest',
                  origin='lower',
                  norm=LogNorm(vmin=max(vmin, 1e-10), vmax=max(vmax, 1e-9)))
        ax.plot(halfwidth, halfwidth, marker='x', mew=5,
                color=vigit_color_1, ms=12)
        ax.set_xlim(0, 2 * halfwidth)
        ax.set_ylim(0, 2 * halfwidth)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        for line in ax.yaxis.get_majorticklines() + ax.xaxis.get_majorticklines():
            line.set_markersize(0)

    if coord_obs[0] is not None:
        ox = float(coord_obs[0]) - cx + halfwidth - 1
        oy = float(coord_obs[1]) - cy + halfwidth - 1
        for ax in (ax_chk, ax_sci):
            ax.errorbar(ox, oy, mew=5, marker='x', color=color_green, ms=12)

    ax_sci.text(0.95, 0.95, '\\textbf{Observed}', ha='right', va='top',
                transform=ax_sci.transAxes, color=color_green,
                fontsize=legend_size,
                path_effects=[PathEffects.withStroke(linewidth=6,
                                                      foreground='w')])
    ax_sci.text(0.95, 0.85, '\\textbf{Expected}', ha='right', va='top',
                transform=ax_sci.transAxes, color=vigit_color_1,
                fontsize=legend_size,
                path_effects=[PathEffects.withStroke(linewidth=6,
                                                      foreground='w')])

    plt.savefig(str(out_dir / (fits_path.stem + '_poststamp.pdf')), dpi=600)
    plt.close()
