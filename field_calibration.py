#!/usr/bin/env python
"""
field_calibration.py — Retrieve and transform photometric reference catalogues.

Queries Vizier (SDSS, PanSTARRS, DES, 2MASS, WISE) and the NOIR DataLab
(Legacy Survey DR10) for a given sky position, applies colour transformations
to produce catalogues in multiple filter systems, and writes them to ASCII
files ready for use with photometry.py.

Usage
-----
    python field_calibration.py --ra 11:33:41.550 --dec 00:43:33.50

Run ``python field_calibration.py --help`` for all options.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import    argparse
from      astroquery.vizier import Vizier
import    astropy.units as u
import    astropy.coordinates as coord
from      astropy.io import ascii
from      astropy import table
import    cat_tools
import    copy
import    fits_tools
from      utils import bcolors
import    numpy as np
from      pathlib import Path
import    sys
import    urllib.request


def get_parser():
    """Build the argument parser for field_calibration.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(
        description=(
            'Retrieve photometric catalogues. '
            '2MASS, PS1, SDSS and SkyMapper are the input catalogs. '
            'Bessel catalogues are generated through colour equations. '
            'PS1, SDSS and SkyMapper are in the AB system; '
            'Bessel and 2MASS are in the Vega system.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--ra',   type=str, required=True,
                   help='RA(J2000) in HMS or decimal degrees. '
                        'For negative Dec write \" -12:20:20.2\"')
    p.add_argument('--dec',  type=str, required=True,
                   help='Dec(J2000) in HMS or decimal degrees.')
    p.add_argument('--radius', type=float, default=10,
                   help='Search radius (arcmin).')
    p.add_argument('--outdir', type=str, default='results/',
                   help='Output directory.')
    p.add_argument('--type',   type=str, default='all',
                   choices=['all', 'optical', 'nir'],
                   help='Catalogue type to generate.')
    return p


def _write_cat(result, outdir, filename, ra_key, dec_key,
               mag_key, err_key, mask_good):
    """Write a filtered, four-column catalogue subset to a headerless ASCII file.

    Parameters
    ----------
    result : astropy.table.Table
        Full catalogue table containing at least the four named columns.
    outdir : pathlib.Path
        Output directory.
    filename : str
        Output filename relative to *outdir*.
    ra_key, dec_key : str
        Column names for right ascension and declination.
    mag_key, err_key : str
        Column names for magnitude and magnitude uncertainty.
    mask_good : array-like of bool
        Boolean mask selecting the rows to include.
    """
    ascii.write(
        result[[ra_key, dec_key, mag_key, err_key]][mask_good],
        str(outdir / filename), overwrite=True, format='no_header')


def main(args=None):
    """Run the field calibration pipeline.

    Parameters
    ----------
    args : list of str, optional
        Command-line arguments (``sys.argv[1:]`` when None).
    """
    parser = get_parser()
    args   = parser.parse_args(args)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ra_dd, dec_dd = fits_tools.convert_hms_dd(
        args.ra, args.dec.replace(' ', '').replace('"', ''))
    coordinates   = coord.SkyCoord(ra_dd, dec_dd, unit=u.deg)

    # ------------------------------------------------------------------
    # SDSS
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\nSDSS catalogues' + bcolors.ENDC)

    try:
        v      = Vizier(row_limit=100000)
        result = v.query_region(
            coordinates, radius=args.radius / 60. * u.deg,
            catalog=cat_tools.catalog_prop['SDSS']['CATID'])
        if len(result) == 0:
            raise ValueError('Field not covered by SDSS')
        result = result[cat_tools.catalog_prop['SDSS']['CATID_OUT']]
        result = result[cat_tools.catalog_prop['SDSS']['KEYWORDS']]
        result = result[(result['class'] == 6) & (result['mode'] == 1)]

        for f in ('u', 'g', 'r', 'i', 'z'):
            result.rename_column(f + 'mag',   f + '_SDSS')
            result.rename_column('e_' + f + 'mag', f + '_SDSS_ERR')

        for key in [k for k in result.colnames if 'SDSS' in k]:
            result[key].format = '.4f'

        result_sdss = copy.deepcopy(result)

        for f in ('u', 'g', 'r', 'i', 'z'):
            mask = ((result[f + '_SDSS_ERR']
                     > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH'])
                    & (result[f + '_SDSS_ERR']
                       < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))
            ascii.write(
                result[['RA_ICRS', 'DE_ICRS',
                         f + '_SDSS', f + '_SDSS_ERR']][mask],
                str(outdir / f'SDSS_SDSS_{f}.cat'),
                overwrite=True, format='no_header')

        # Bessel
        result = cat_tools.sdss_to_bessel(result_sdss)
        for key in [k for k in result.colnames if 'BESSEL' in k]:
            result[key].format = '.4f'
        for f in ('B', 'V', 'R', 'I'):
            mask = ((result[f + '_BESSEL_ERR'] > 0.)
                    & (result[f + '_BESSEL_ERR'] < 0.3))
            ascii.write(
                result[['RA_ICRS', 'DE_ICRS',
                         f + '_BESSEL', f + '_BESSEL_ERR']][mask],
                str(outdir / f'SDSS_BESSEL_{f}.cat'),
                overwrite=True, format='no_header')

        # GROND
        print(bcolors.WARNING + 'Convert SDSS -> GROND' + bcolors.ENDC)
        result = cat_tools.sdss_to_grond(result_sdss)
        for key in [k for k in result.colnames if 'GROND' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r', 'i', 'z'):
            mask = ((result[f + '_GROND_ERR'] > 0.)
                    & (result[f + '_GROND_ERR'] < 0.3))
            ascii.write(
                result[['RA_ICRS', 'DE_ICRS',
                         f + '_GROND', f + '_GROND_ERR']][mask],
                str(outdir / f'SDSS_GROND_{f}.cat'),
                overwrite=True, format='no_header')

        # DES  (was printing "Convert SDSS -> GROND" — fixed)
        print(bcolors.WARNING + 'Convert SDSS -> DES' + bcolors.ENDC)
        result = cat_tools.sdss_to_des(result_sdss)
        for key in [k for k in result.colnames if 'DES' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r', 'i', 'z'):
            mask = ((result[f + '_DES_ERR'] > 0.)
                    & (result[f + '_DES_ERR'] < 0.3))
            ascii.write(
                result[['RA_ICRS', 'DE_ICRS',
                         f + '_DES', f + '_DES_ERR']][mask],
                str(outdir / f'SDSS_DES_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.OKGREEN + 'SDSS query successful.' + bcolors.ENDC)

    except Exception as exc:
        print(bcolors.FAIL
              + f'Field not covered by SDSS or Vizier query failed: {exc}'
              + bcolors.ENDC)

    # ------------------------------------------------------------------
    # PanSTARRS (cross-matched with Gaia for star/galaxy separation)
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\nPS1 catalogues' + bcolors.ENDC)

    try:
        v           = Vizier(columns=['all'], row_limit=-1)
        result_gaia = v.query_region(
            coordinates, radius=args.radius * u.arcmin,
            catalog=cat_tools.catalog_prop['GAIA']['CATID'])
        result_gaia = result_gaia[cat_tools.catalog_prop['GAIA']['CATID_OUT']]
        result_gaia = result_gaia[cat_tools.catalog_prop['GAIA']['KEYWORDS']]
        result_gaia = result_gaia[
            (result_gaia['Plx'] > 0)
            & (result_gaia['e_pmRA'] > 0)
            & (result_gaia['e_pmDE'] > 0)]
        result_gaia = result_gaia[
            (result_gaia['Plx'] / result_gaia['e_Plx'] > 3)
            | (abs(result_gaia['pmRA'] / result_gaia['e_pmRA']) > 3)
            | (abs(result_gaia['pmDE'] / result_gaia['e_pmDE']) > 3)]

        v              = Vizier(columns=['all'], row_limit=-1)
        result_ps1     = v.query_region(
            coordinates, radius=args.radius * u.arcmin,
            catalog=cat_tools.catalog_prop['PanSTARRS']['CATID'])
        result_ps1     = result_ps1[cat_tools.catalog_prop['PanSTARRS']['CATID_OUT']]
        ps1_keys       = ['RAJ2000', 'DEJ2000',
                          'gmag', 'e_gmag', 'rmag', 'e_rmag',
                          'imag', 'e_imag', 'zmag', 'e_zmag',
                          'ymag', 'e_ymag']
        result_ps1     = result_ps1[ps1_keys]

        # Cross-match Gaia × PS1 to keep only confirmed stars
        gaia_arr = np.array(list(np.asarray(result_gaia).tolist()))
        ps1_arr  = np.array(list(np.asarray(result_ps1).tolist()))

        matched_arr  = cat_tools.crossmatch_catalogs(gaia_arr, ps1_arr, 0.5)
        matched      = table.Table(
            matched_arr,
            names=ps1_keys
            + ['GAIA_' + k for k in cat_tools.catalog_prop['GAIA']['KEYWORDS']]
            + ['DIST'])
        matched['DIST'] *= 3600.
        matched          = matched[ps1_keys]

        for f in ('g', 'r', 'i', 'z', 'y'):
            matched.rename_column(f + 'mag',   f + '_PS1')
            matched.rename_column('e_' + f + 'mag', f + '_PS1_ERR')

        for key in [k for k in matched.colnames if 'PS1' in k]:
            matched[key].format = '.4f'

        result_ps1_clean = copy.deepcopy(matched)

        # Write PS1
        for f in ('g', 'r', 'i', 'z', 'y'):
            mask = ((matched[f + '_PS1_ERR']
                     > cat_tools.catalog_prop['PanSTARRS']['SIGMA_HIGH'])
                    & (matched[f + '_PS1_ERR']
                       < cat_tools.catalog_prop['PanSTARRS']['SIGMA_LOW']))
            ascii.write(
                matched[['RAJ2000', 'DEJ2000',
                          f + '_PS1', f + '_PS1_ERR']][mask],
                str(outdir / f'PS1_PS1_{f}.cat'),
                overwrite=True, format='no_header')

        # PS1 → SDSS
        print(bcolors.WARNING + 'Convert PS1 -> SDSS' + bcolors.ENDC)
        result        = cat_tools.ps1_to_sdss(result_ps1_clean)
        for key in [k for k in result.colnames if 'SDSS' in k]:
            result[key].format = '.4f'
        result_sdss_from_ps1 = copy.deepcopy(result)
        for f in ('u', 'g', 'r', 'i', 'z'):
            mask = ((result[f + '_SDSS_ERR']
                     > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH'])
                    & (result[f + '_SDSS_ERR']
                       < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_SDSS', f + '_SDSS_ERR']][mask],
                str(outdir / f'PS1_SDSS_{f}.cat'),
                overwrite=True, format='no_header')

        # PS1 → HSC
        print(bcolors.WARNING + 'Convert PS1 -> HSC' + bcolors.ENDC)
        result = cat_tools.ps1_to_hsc(result_ps1_clean)
        for key in [k for k in result.colnames if 'HSC' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r', 'i', 'z', 'y'):
            mask = (result[f + '_HSC_ERR'] > 0.) & (result[f + '_HSC_ERR'] < 0.3)
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_HSC', f + '_HSC_ERR']][mask],
                str(outdir / f'PS1_HSC_{f}.cat'),
                overwrite=True, format='no_header')

        # PS1 → SDSS → Bessel
        print(bcolors.WARNING + 'Convert PS1 -> SDSS -> Bessel' + bcolors.ENDC)
        result = cat_tools.sdss_to_bessel(result_sdss_from_ps1)
        for key in [k for k in result.colnames if 'BESSEL' in k]:
            result[key].format = '.4f'
        for f in ('B', 'V', 'R', 'I'):
            mask = (result[f + '_BESSEL_ERR'] > 0.) & (result[f + '_BESSEL_ERR'] < 0.3)
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_BESSEL', f + '_BESSEL_ERR']][mask],
                str(outdir / f'PS1_BESSEL_{f}.cat'),
                overwrite=True, format='no_header')

        # PS1 → ZTF
        print(bcolors.WARNING + 'Convert PS1 -> ZTF' + bcolors.ENDC)
        result = cat_tools.ps1_to_ztf(result_ps1_clean)
        for key in [k for k in result.colnames if 'ZTF' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r'):
            mask = (result[f + '_ZTF_ERR'] > 0.) & (result[f + '_ZTF_ERR'] < 0.3)
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_ZTF', f + '_ZTF_ERR']][mask],
                str(outdir / f'PS1_ZTF_{f}.cat'),
                overwrite=True, format='no_header')

        # PS1 → SDSS → GROND
        print(bcolors.WARNING + 'Convert PS1 -> SDSS -> GROND' + bcolors.ENDC)
        result = cat_tools.sdss_to_grond(result_sdss_from_ps1)
        for key in [k for k in result.colnames if 'GROND' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r', 'i', 'z'):
            mask = (result[f + '_GROND_ERR'] > 0.) & (result[f + '_GROND_ERR'] < 0.3)
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_GROND', f + '_GROND_ERR']][mask],
                str(outdir / f'PS1_GROND_{f}.cat'),
                overwrite=True, format='no_header')

        # PS1 → SDSS → DES
        print(bcolors.WARNING + 'Convert PS1 -> SDSS -> DES' + bcolors.ENDC)
        result = cat_tools.sdss_to_des(result_sdss_from_ps1)
        for key in [k for k in result.colnames if 'DES' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r', 'i', 'z'):
            mask = (result[f + '_DES_ERR'] > 0.) & (result[f + '_DES_ERR'] < 0.3)
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_DES', f + '_DES_ERR']][mask],
                str(outdir / f'PS1_DES_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.OKGREEN + 'PS1 query successful.' + bcolors.ENDC)

    except Exception as exc:
        print(bcolors.FAIL
              + f'Field not covered by PS1 or query failed: {exc}'
              + bcolors.ENDC)

    # ------------------------------------------------------------------
    # DESI Legacy Imaging Survey DR10 (via NOIR DataLab)
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\nDESI Legacy Imaging Surveys DR10' + bcolors.ENDC)

    try:
        import routines_noir
        result = cat_tools.query_noir_datalab_stars(
            coordinates.ra.value, coordinates.dec.value, args.radius)

        if len(result) > 0:
            for f in ('g', 'r', 'i', 'z'):
                temp = result['ra', 'dec', 'mag_' + f, 'snr_' + f,
                               'allmask_' + f].copy()
                temp['ra'].format  = '.6f'
                temp['dec'].format = '.6f'
                mask = (temp['snr_' + f] >= 3) & (temp['allmask_' + f] == 0)
                temp = temp[mask]
                if len(temp):
                    temp['err_' + f] = 2.5 * np.log10(1. + 1. / temp['snr_' + f])
                    for k in ('mag_' + f, 'err_' + f):
                        temp[k].format = '.4f'
                    ascii.write(temp['ra', 'dec', 'mag_' + f, 'err_' + f],
                                str(outdir / f'LS_LS_{f}.cat'),
                                overwrite=True, format='no_header')

            print(bcolors.WARNING + 'LS query successful.' + bcolors.ENDC)

            # LS → SDSS → Bessel
            qual_mask = (
                (result['allmask_g'] == 0) & (result['snr_g'] >= 3)
                & (result['allmask_r'] == 0) & (result['snr_r'] >= 3)
                & (result['allmask_z'] == 0) & (result['snr_z'] >= 3))
            result = result[qual_mask]

            for f in ('g', 'r', 'z'):
                result['err_' + f] = 2.5 * np.log10(1. + 1. / result['snr_' + f])
                result.rename_column('mag_' + f, f + '_DES')
                result.rename_column('err_' + f, f + '_DES_ERR')

            for key in [k for k in result.colnames if 'DES' in k]:
                result[key].format = '.4f'

            result['RAJ2000'] = result['ra']
            result['DEJ2000'] = result['dec']
            result['i_DES']     = result['z_DES']
            result['i_DES_ERR'] = result['z_DES_ERR']

            print(bcolors.WARNING + 'Convert LS -> SDSS' + bcolors.ENDC)
            result = cat_tools.des_to_sdss(result)
            for key in [k for k in result.colnames if 'SDSS' in k]:
                result[key].format = '.4f'
            result_sdss_from_ls = copy.deepcopy(result)

            for f in ('g', 'r', 'i'):
                mask = ((result[f + '_SDSS_ERR']
                         > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH'])
                        & (result[f + '_SDSS_ERR']
                           < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))
                ascii.write(
                    result[['RAJ2000', 'DEJ2000',
                             f + '_SDSS', f + '_SDSS_ERR']][mask],
                    str(outdir / f'LS_SDSS_{f}.cat'),
                    overwrite=True, format='no_header')

            print(bcolors.WARNING + 'Convert LS -> SDSS -> Bessel' + bcolors.ENDC)
            result = cat_tools.sdss_to_bessel(result_sdss_from_ls)
            for key in [k for k in result.colnames if 'BESSEL' in k]:
                result[key].format = '.4f'
            for f in ('B', 'V', 'R', 'I'):
                mask = ((result[f + '_BESSEL_ERR'] > 0.)
                        & (result[f + '_BESSEL_ERR'] < 0.3))
                ascii.write(
                    result[['RAJ2000', 'DEJ2000',
                             f + '_BESSEL', f + '_BESSEL_ERR']][mask],
                    str(outdir / f'LS_BESSEL_{f}.cat'),
                    overwrite=True, format='no_header')
        else:
            print(bcolors.FAIL + 'Field not covered by Legacy Survey.' + bcolors.ENDC)

    except ImportError:
        print(bcolors.WARNING
              + 'routines_noir not available (dl library not installed). '
              'Skipping Legacy Survey query.'
              + bcolors.ENDC)
    except Exception as exc:
        print(bcolors.FAIL + f'NOIR DataLab query failed: {exc}' + bcolors.ENDC)

    # ------------------------------------------------------------------
    # DES (via Vizier)
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\nDES catalogues' + bcolors.ENDC)

    try:
        v      = Vizier(columns=['all'], row_limit=-1)
        result = v.query_region(
            coordinates, radius=args.radius / 60. * u.deg,
            catalog=cat_tools.catalog_prop['DES']['CATID'])
        if len(result) == 0:
            raise ValueError('Field not covered by DES')
        result = result[cat_tools.catalog_prop['DES']['CATID_OUT']]
        result = result[cat_tools.catalog_prop['DES']['KEYWORDS']]
        result = result[
            (result['S_Gg'] > 0.85) & (result['S_Gr'] > 0.85)
            & (result['S_Gi'] > 0.85)
            & (result['gFlag'] <= 4) & (result['rFlag'] <= 4)
            & (result['iFlag'] <= 4) & (result['zFlag'] <= 4)]

        for f in ('g', 'r', 'i', 'z', 'Y'):
            result.rename_column(f + 'mag',   f + '_DES')
            result.rename_column('e_' + f + 'mag', f + '_DES_ERR')

        for key in [k for k in result.colnames if 'DES' in k]:
            result[key].format = '.4f'

        des_raw = copy.deepcopy(result)

        for f in ('g', 'r', 'i', 'z', 'Y'):
            mask = ((result[f + '_DES_ERR']
                     > cat_tools.catalog_prop['DES']['SIGMA_HIGH'])
                    & (result[f + '_DES_ERR']
                       < cat_tools.catalog_prop['DES']['SIGMA_LOW']))
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_DES', f + '_DES_ERR']][mask],
                str(outdir / f'DES_DES_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.WARNING + 'Convert DES -> SDSS' + bcolors.ENDC)
        result = cat_tools.des_to_sdss(des_raw)
        for key in [k for k in result.colnames if 'SDSS' in k]:
            result[key].format = '.4f'
        result_sdss_from_des = copy.deepcopy(result)
        for f in ('g', 'r', 'i', 'z'):
            mask = ((result[f + '_SDSS_ERR']
                     > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH'])
                    & (result[f + '_SDSS_ERR']
                       < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_SDSS', f + '_SDSS_ERR']][mask],
                str(outdir / f'DES_SDSS_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.WARNING + 'Convert DES -> SDSS -> Bessel' + bcolors.ENDC)
        result = cat_tools.sdss_to_bessel(result_sdss_from_des)
        for key in [k for k in result.colnames if 'BESSEL' in k]:
            result[key].format = '.4f'
        for f in ('B', 'V', 'R', 'I'):
            mask = ((result[f + '_BESSEL_ERR'] > 0.)
                    & (result[f + '_BESSEL_ERR'] < 0.3))
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_BESSEL', f + '_BESSEL_ERR']][mask],
                str(outdir / f'DES_BESSEL_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.WARNING + 'Convert DES -> SDSS -> GROND' + bcolors.ENDC)
        result = cat_tools.sdss_to_grond(result_sdss_from_des)
        for key in [k for k in result.colnames if 'GROND' in k]:
            result[key].format = '.4f'
        for f in ('g', 'r', 'i', 'z'):
            mask = ((result[f + '_GROND_ERR'] > 0.)
                    & (result[f + '_GROND_ERR'] < 0.3))
            ascii.write(
                result[['RAJ2000', 'DEJ2000',
                         f + '_GROND', f + '_GROND_ERR']][mask],
                str(outdir / f'DES_GROND_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.OKGREEN + 'DES query successful.' + bcolors.ENDC)

    except Exception as exc:
        print(bcolors.FAIL + f'Field not covered by DES: {exc}' + bcolors.ENDC)

    # ------------------------------------------------------------------
    # SkyMapper (southern sky only)
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\nSkyMapper catalogue' + bcolors.ENDC)

    if coordinates.dec.deg < 0:
        try:
            url = (
                'http://skymapper.anu.edu.au/sm-cone/public/query'
                f'?RA={coordinates.ra.deg:.3f}'
                f'&DEC={coordinates.dec.deg:.3f}'
                f'&SR={args.radius / 60.:.1f}'
                '&RESPONSEFORMAT=CSV')
            raw = urllib.request.urlopen(url).read()
            csv_path = outdir / 'skymapper.csv'
            csv_path.write_bytes(raw)

            cat_sm = ascii.read(str(csv_path), format='csv')
            csv_path.unlink()

            cat_sm = cat_sm[
                (cat_sm['class_star'] > 0.95) & (cat_sm['flags'] == 0)
                & (cat_sm['g_psf'] > 0) & (cat_sm['r_psf'] > 0)
                & (cat_sm['i_psf'] > 0) & (cat_sm['z_psf'] > 0)]

            sm_keys = ['raj2000', 'dej2000',
                       'u_psf', 'e_u_psf', 'v_psf', 'e_v_psf',
                       'g_psf', 'e_g_psf', 'r_psf', 'e_r_psf',
                       'i_psf', 'e_i_psf', 'z_psf', 'e_z_psf']
            cat_sm = cat_sm[sm_keys]
            cat_sm.rename_column('raj2000', 'RAJ2000')
            cat_sm.rename_column('dej2000', 'DEJ2000')

            for f in ('u', 'g', 'r', 'i', 'z'):
                cat_sm.rename_column(f + '_psf',   f + '_SDSS')
                cat_sm.rename_column('e_' + f + '_psf', f + '_SDSS_ERR')

            for key in [k for k in cat_sm.colnames if 'SDSS' in k]:
                cat_sm[key].format = '.4f'

            for f in ('u', 'g', 'r', 'i', 'z'):
                mask = ((cat_sm[f + '_SDSS_ERR']
                         > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH'])
                        & (cat_sm[f + '_SDSS_ERR']
                           < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))
                ascii.write(
                    cat_sm[['RAJ2000', 'DEJ2000',
                             f + '_SDSS', f + '_SDSS_ERR']][mask],
                    str(outdir / f'SkyMapper_SDSS_{f}.cat'),
                    overwrite=True, format='no_header')

            cat_sm = cat_tools.sdss_to_bessel(cat_sm)
            for key in [k for k in cat_sm.colnames if 'BESSEL' in k]:
                cat_sm[key].format = '.4f'
            for f in ('B', 'V', 'R', 'I'):
                mask = ((cat_sm[f + '_BESSEL_ERR'] > 0.)
                        & (cat_sm[f + '_BESSEL_ERR'] < 0.3))
                ascii.write(
                    cat_sm[['RAJ2000', 'DEJ2000',
                             f + '_BESSEL', f + '_BESSEL_ERR']][mask],
                    str(outdir / f'SkyMapper_BESSEL_{f}.cat'),
                    overwrite=True, format='no_header')

            print(bcolors.OKGREEN + 'SkyMapper query successful.' + bcolors.ENDC)

        except Exception as exc:
            print(bcolors.WARNING
                  + f'SkyMapper query failed: {exc}'
                  + bcolors.ENDC)
    else:
        print(bcolors.FAIL
              + 'Declination is above 0°. No SkyMapper data.'
              + bcolors.ENDC)

    # ------------------------------------------------------------------
    # 2MASS
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\n2MASS catalogues' + bcolors.ENDC)

    try:
        v      = Vizier(row_limit=10000)
        result = v.query_region(
            coord.SkyCoord(ra_dd, dec_dd, unit=u.deg),
            radius=args.radius / 60. * u.deg,
            catalog=cat_tools.catalog_prop['2MASS']['CATID'])
        if len(result) == 0:
            raise ValueError('Field not covered by 2MASS')
        result = result[cat_tools.catalog_prop['2MASS']['CATID_OUT']]
        result = result[cat_tools.catalog_prop['2MASS']['KEYWORDS']]

        for key in [k for k in result.colnames if 'mag' in k]:
            result[key].format = '.4f'

        for f in ('J', 'H', 'K'):
            mask = ((result['e_' + f + 'mag']
                     > cat_tools.catalog_prop['2MASS']['SIGMA_HIGH'])
                    & (result['e_' + f + 'mag']
                       < cat_tools.catalog_prop['2MASS']['SIGMA_LOW']))
            ascii.write(
                result[['RAJ2000', 'DEJ2000', f + 'mag', 'e_' + f + 'mag']][mask],
                str(outdir / f'2MASS_{f}.cat'),
                overwrite=True, format='no_header')

        print(bcolors.OKGREEN + '2MASS query successful.' + bcolors.ENDC)

    except Exception as exc:
        print(bcolors.FAIL + f'2MASS Vizier query failed: {exc}' + bcolors.ENDC)

    # ------------------------------------------------------------------
    # WISE (AllWISE)
    # ------------------------------------------------------------------
    print(bcolors.OKBLUE + '\nAllWISE' + bcolors.ENDC)

    try:
        v         = Vizier(columns=['all'],
                           catalog=cat_tools.catalog_prop['WISE']['CATID'],
                           row_limit=-1)
        data_wise = v.query_region(coordinates, radius=50 * u.arcmin)
        data_wise = data_wise[cat_tools.catalog_prop['WISE']['CATID_OUT']]
        data_wise = data_wise[cat_tools.catalog_prop['WISE']['KEYWORDS']]

        mask = ((data_wise['e_W1mag'] < 0.3) & (data_wise['e_W1mag'] > 0.)
                & (data_wise['e_W2mag'] < 0.3) & (data_wise['e_W2mag'] > 0.))
        data_wise = data_wise[mask]

        ascii.write(data_wise[['RAJ2000', 'DEJ2000', 'W1mag', 'e_W1mag']],
                    str(outdir / 'unWISE_W1.cat'),
                    format='no_header', overwrite=True)
        ascii.write(data_wise[['RAJ2000', 'DEJ2000', 'W2mag', 'e_W2mag']],
                    str(outdir / 'unWISE_W2.cat'),
                    format='no_header', overwrite=True)
        ascii.write(data_wise, str(outdir / 'unWISE.cat'), overwrite=True)

        print(bcolors.OKGREEN + 'WISE query successful.' + bcolors.ENDC)

    except Exception as exc:
        print(bcolors.FAIL + f'WISE Vizier query failed: {exc}' + bcolors.ENDC)

    print('\n')


if __name__ == '__main__':
    main(sys.argv[1:])
