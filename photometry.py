#!/usr/bin/env python
"""
photometry.py — Aperture photometry pipeline for ground-based imaging.

Workflow
--------
1. Verify that the target lies in the image footprint (WCS check).
2. Retrieve or read a reference photometric catalogue.
3. Build a local comparison-star sequence (sep + cross-match).
4. Determine the photometric zeropoint (bootstrap, per aperture).
5. Run aperture photometry on the science target.
6. Output calibrated magnitudes and diagnostic plots.

Usage
-----
    python photometry.py --ra 11:33:41.550 --dec 00:43:33.50 \\
        --fits IMAGE.fits --ref-file results/SDSS_SDSS_r.cat

Run ``python photometry.py --help`` for all options.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import   argparse
import   copy
import   glob
import   logging
import   shutil
import   sys
from     pathlib import Path

import   numpy as np
import   numpy.lib.recfunctions as rfn
from     astropy import table, units as u
from     astropy.io import ascii, fits
from     matplotlib import pylab as plt

import   calibration
import   cat_tools
import   extraction
import   fits_tools
from     plotsettings import *
from     utils import bcolors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    p = argparse.ArgumentParser(
        description='Aperture photometry pipeline for ground-based images.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--ra',  type=str, required=True,
                   help='RA(J2000) in HMS or decimal degrees.')
    p.add_argument('--dec', type=str, required=True,
                   help='Dec(J2000) in HMS or decimal degrees.')
    p.add_argument('--fits', type=str, required=True,
                   help='Science FITS filename.')
    p.add_argument('--host-offset', type=float, default=5,
                   help='Detection search radius around target (arcsec).')

    p.add_argument('--ref-cat',    type=str, default=None,
                   help='Reference catalogue name (SDSS, 2MASS, PS1, …).')
    p.add_argument('--ref-filter', type=str, default=None,
                   help='Filter of the reference catalogue.')
    p.add_argument('--ref-file',   type=str, default=None,
                   help='Local reference catalogue file. Overrides --ref-cat.')
    p.add_argument('--ref-image',  type=str, default='',
                   help='Reference image for sep dual-image mode.')
    p.add_argument('--ref-radius', type=float, default=10,
                   help='Vizier search radius (arcmin).')

    p.add_argument('--ana-thresh',       type=float, default=1,
                   help='Analysis threshold (sigma).')
    p.add_argument('--ap-diam',          type=float, default=[1., 1.5, 2., 3.],
                   nargs='+', help='Aperture diameters as multiples of FWHM.')
    p.add_argument('--ap-diam-ul',       type=float, default=2, nargs='+',
                   help='Aperture diameter for forced-phot upper limit.')
    p.add_argument('--det-thresh',       type=float, default=1,
                   help='Detection threshold (sigma).')
    p.add_argument('--gain',             type=str,   default=None,
                   help='Header keyword for CCD gain (e-/ADU).')
    p.add_argument('--back-size',        type=int,   default=64)
    p.add_argument('--back-filtersize',  type=int,   default=3)
    p.add_argument('--deblend-nthresh',  type=int,   default=64)
    p.add_argument('--deblend-mincont',  type=float, default=0.0001)
    p.add_argument('--noflags',          action='store_true', default=False)

    p.add_argument('--mag-cut',       type=float, default=12.)
    p.add_argument('--mag-stdfaint',  type=str,   default='0')
    p.add_argument('--mag-stdbright', type=str,   default='0')
    p.add_argument('--maxstars',      type=int,   default=200)

    p.add_argument('--auto',      action='store_true', default=False)
    p.add_argument('--bw',        action='store_true', default=False,
                   help='Monochrome terminal output.')
    p.add_argument('--keeptemp',  action='store_true', default=False)
    p.add_argument('--loglevel',  type=str, default='INFO')
    p.add_argument('--outdir',    type=str, default='results/')
    p.add_argument('--sex-loglevel', type=str, default='WARNING')
    p.add_argument('--tol',       type=float, default=1,
                   help='Cross-matching tolerance (arcsec).')

    return p


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    """Run the aperture photometry pipeline."""
    parser = get_parser()
    args   = parser.parse_args(argv)

    args.mag_stdfaint  = float(args.mag_stdfaint)
    args.mag_stdbright = float(args.mag_stdbright)

    # Logger
    log_path = Path(args.fits).with_suffix('.log')
    logger   = logging.getLogger('photometry')
    logger.setLevel(args.loglevel)
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(args.loglevel)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    # -----------------------------------------------------------------------
    # Step 1: Administration
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 1: Administration\n' + bcolors.ENDC)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info(f'Output directory: {outdir}')

    ra_dd, dec_dd = fits_tools.convert_hms_dd(args.ra, args.dec)

    object_properties = table.Table()
    object_properties['OBJECT']  = [Path(args.fits).stem]
    object_properties['RA_HMS']  = [args.ra]
    object_properties['DEC_DMS'] = [args.dec]
    object_properties['RA']      = [ra_dd]
    object_properties['DEC']     = [dec_dd]

    print(bcolors.OKGREEN + '\nIs the target in the image footprint?' + bcolors.ENDC)
    result = fits_tools.sky2xy(args.fits, ra=ra_dd, dec=dec_dd)
    if result is False:
        msg = (f'Target (RA={args.ra}, Dec={args.dec}) is outside the '
               f'footprint of {args.fits}.')
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    x_exp, y_exp = result
    object_properties['X_EXP'] = [x_exp]
    object_properties['Y_EXP'] = [y_exp]

    science_xy = outdir / (object_properties['OBJECT'][0] + '_xy.cat')
    ascii.write(np.array([[x_exp, y_exp]]),
                str(science_xy), format='no_header', overwrite=True)
    print('Yes.')
    logger.info('Target in footprint: OK')

    pix_scale = fits_tools.pix2arcsec(args.fits)

    # -----------------------------------------------------------------------
    # Step 2: Reference catalogue
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 2: Flux calibration\n' + bcolors.ENDC)

    if args.ref_file is not None:
        msg = f'Using catalogue: {args.ref_file}'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        ref_cat = ascii.read(args.ref_file, names=('RA', 'DEC', 'MAG', 'MAGERR'))
        object_properties['PHOTCAL'] = [args.ref_file]

        filename_stars = outdir / (object_properties['OBJECT'][0]
                                   + '_' + Path(args.ref_file).stem + '.cat')
        shutil.copy(args.ref_file,
                    str(filename_stars).replace('.cat', '_refcat.cat'))

    else:
        msg = 'Retrieving catalogue from Vizier'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        object_properties['PHOTCAL'] = [args.ref_cat + '/' + args.ref_filter]

        filename_stars = outdir / (object_properties['OBJECT'][0]
                                   + '_' + args.ref_cat
                                   + '_' + args.ref_filter + '.cat')
        cat_tools.retrieve_photcat(
            object_properties, args.ref_cat, cat_tools.catalog_prop,
            filename=str(filename_stars),
            radius=args.ref_radius * u.arcmin,
            output_dir=str(outdir))

        shutil.copy(str(filename_stars),
                    str(filename_stars).replace('.cat', '_refcat.cat'))

    # Build local sequence
    print(bcolors.OKGREEN + bcolors.BOLD + '\nBuilding the local sequence' + bcolors.ENDC)

    refcat = str(filename_stars).replace('.cat', '_refcat.cat')
    ref_cat = ascii.read(refcat, names=('RA', 'DEC', 'MAG', 'MAGERR'))

    ref_cat = ref_cat[ref_cat['MAG'] > args.mag_cut]
    logger.info(f'Removing stars brighter than {args.mag_cut:.2f} mag')
    ascii.write(ref_cat, refcat, overwrite=True, format='no_header')

    refcat_xy = refcat.replace('_refcat.cat', '_refcat_xy.cat')
    catalog_clean = fits_tools.sky2xy(args.fits, catalog=refcat)

    if not len(catalog_clean):
        msg = 'No reference stars in the image footprint.'
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    ascii.write(catalog_clean, refcat_xy, overwrite=True, format='no_header')
    print(bcolors.OKGREEN
          + f'\n{len(catalog_clean)} reference stars in footprint.' + bcolors.ENDC)

    # -----------------------------------------------------------------------
    # Run sep on reference stars (local sequence detection)
    # -----------------------------------------------------------------------
    print(bcolors.OKGREEN + '\nSelecting stars for local sequence' + bcolors.ENDC)

    ref_stars = extraction.extract_sources(
        args.fits,
        analysis_thresh  = 3,
        assoc_file       = refcat_xy,
        assoc_cols       = "1,2",
        assoc_radius     = args.tol / pix_scale,
        detect_thresh    = 3,
        flag             = 'loc_seq',
        gain_key         = args.gain,
        logger           = logger,
        output_dir       = str(outdir),
        phot_apertures   = np.array([10.]),
        ref_image        = args.ref_image)

    if not len(ref_stars):
        msg = ('No reference stars detected. Check --tol, image WCS, '
               'and catalogue.')
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    ref_stars = ref_stars[
        (ref_stars['MAGERR_APER'] > 0.0002)
        & (ref_stars['MAGERR_AUTO'] > 0.0002)]
    logger.info(f'{len(ref_stars)} stars for local sequence')

    def _prune_stars(stars: table.Table) -> table.Table:
        if not args.noflags:
            stars = stars[stars['FLAGS'] == 0]
        fwhm = np.asarray(stars['FWHM_IMAGE'])
        if np.any(fwhm > 0):
            p25, p75 = np.percentile(fwhm[fwhm > 0], [25, 75])
            stars    = stars[fwhm < p75 + 1.5 * (p75 - p25)]
        if len(stars) > args.maxstars:
            print(bcolors.WARNING
                  + f'Truncating to brightest {args.maxstars} stars.'
                  + bcolors.ENDC)
            stars.sort('MAGERR_AUTO')
            stars = stars[:args.maxstars]
        return stars

    def _crossmatch_to_refcat(stars: table.Table) -> table.Table:
        keys_s = [k for k in (
            'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE',
            'MAG_AUTO', 'MAGERR_AUTO', 'MAG_PETRO', 'MAGERR_PETRO',
            'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_PETRO', 'FLUXERR_PETRO',
            'FWHM_IMAGE', 'FWHM_WORLD', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
            'KRON_RADIUS', 'FLUX_RADIUS', 'FLAGS', 'VECTOR_ASSOC',
            'NUMBER_ASSOC', 'MAG_APER', 'MAGERR_APER', 'FLUX_APER',
            'FLUXERR_APER') if k in stars.colnames]
        keys_r = list(ref_cat.colnames)

        arr_s = rfn.structured_to_unstructured(
            np.asarray(stars[keys_s]), dtype=float)
        arr_r = rfn.structured_to_unstructured(
            np.asarray(ref_cat[keys_r]), dtype=float)

        matched = cat_tools.crossmatch_catalogs(arr_s, arr_r, args.tol)
        if not len(matched):
            msg = (f'No cross-matches within {args.tol} arcsec. '
                   'Check --tol and WCS.')
            print(bcolors.FAIL + msg + bcolors.ENDC)
            logger.error(msg)
            sys.exit(1)

        tbl                   = table.Table(matched, names=keys_r + keys_s + ['DIST'])
        tbl['DIST']           = tbl['DIST'] * 3600.0
        tbl['DIST'].format    = '.3f'
        tbl.rename_column('MAG_APER',    'MAG_INS')
        tbl.rename_column('MAGERR_APER', 'MAGERR_INS')
        tbl.rename_column('MAG',         'MAG_CAT')
        tbl.rename_column('MAGERR',      'MAGERR_CAT')
        tbl['MAG_ZP_TEMP']        = tbl['MAG_CAT'] - tbl['MAG_INS']
        tbl['MAG_ZP_TEMP'].format = '.3f'
        return tbl

    if args.mag_stdbright == 0 and args.mag_stdfaint == 0:
        ref_stars    = _prune_stars(ref_stars)
        logger.info(f'{len(ref_stars)} stars remain after pruning')

        msg = 'Cross-matching catalogues'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        matched_standard = _crossmatch_to_refcat(ref_stars)

        temp = copy.deepcopy(matched_standard)
        temp.sort('MAG_CAT')
        print(temp[['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE',
                     'YWIN_IMAGE', 'FWHM_IMAGE', 'MAG_CAT', 'MAGERR_CAT',
                     'MAG_INS', 'MAGERR_INS', 'MAG_ZP_TEMP', 'DIST']])

        local_seq = calibration.local_sequence(
            matched_standard,
            auto        = args.auto,
            output_file = str(filename_stars).replace('.cat', '_loc_cleaned.cat'),
            fits_path   = args.fits,
            logger      = logger,
            lower       = 5,
            output_dir  = str(outdir),
            upper       = 90)

    else:
        ref_stars = ref_stars[
            (ref_stars['MAG_APER'] >= args.mag_stdbright)
            & (ref_stars['MAG_APER'] <= args.mag_stdfaint)]
        ref_stars = _prune_stars(ref_stars)

        msg = 'Cross-matching catalogues'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        matched_standard = _crossmatch_to_refcat(ref_stars)

        local_seq = calibration.local_sequence(
            matched_standard,
            auto        = args.auto,
            output_file = str(filename_stars).replace('.cat', '_loc_cleaned.cat'),
            fits_path   = args.fits,
            logger      = logger,
            lower       = 0,
            output_dir  = str(outdir),
            upper       = 100)

    matched_standard = local_seq['CAT']
    logger.info(f'Local sequence: {local_seq["NUMSTARS"]} stars')
    print(bcolors.OKGREEN
          + f'\nFinal local sequence: {local_seq["NUMSTARS"]} stars'
          + bcolors.ENDC)

    # FWHM for aperture scaling
    fwhm_vals = np.asarray(matched_standard['FWHM_IMAGE'])
    if np.all(fwhm_vals == 0):
        print(bcolors.WARNING + 'All FWHM_IMAGE = 0; using FLUX_RADIUS.' + bcolors.ENDC)
        FWHM_median = float(np.median(matched_standard['FLUX_RADIUS'])) * 2.0 / 1.1
    else:
        ok          = fwhm_vals[fwhm_vals > 0]
        FWHM_median = float(np.median(ok)) if len(ok) else 3.0

    apertures = 2.0 * FWHM_median * np.array(args.ap_diam)
    logger.info(f'FWHM = {FWHM_median:.2f} px; apertures = {apertures} px')

    # -----------------------------------------------------------------------
    # Step 3: Zeropoint
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 3: Zeropoint calculation\n' + bcolors.ENDC)

    phot_stars = extraction.extract_sources(
        args.fits,
        analysis_thresh  = 3,
        assoc_file       = str(filename_stars).replace('.cat', '_loc_cleaned.cat'),
        assoc_cols       = "1,2",
        assoc_radius     = args.tol / pix_scale,
        detect_thresh    = 3,
        flag             = 'ref_star',
        gain_key         = args.gain,
        logger           = logger,
        output_dir       = str(outdir),
        phot_apertures   = apertures,
        ref_image        = args.ref_image)

    summary_zeropoint = calibration.zeropoint(
        matched_standard, phot_stars,
        fits_path  = args.fits,
        logger     = logger,
        niter      = 1000,
        output_dir = str(outdir),
        tolerance  = args.tol)

    summary_zeropoint['r(FWHM)'] = np.nan

    # -----------------------------------------------------------------------
    # Step 4: Science photometry
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 4: Aperture photometry\n' + bcolors.ENDC)

    phot_all = extraction.extract_sources(
        args.fits,
        analysis_thresh  = 1,
        back_size        = args.back_size,
        back_filtersize  = args.back_filtersize,
        deblend_nthresh  = args.deblend_nthresh,
        deblend_mincont  = args.deblend_mincont,
        detect_thresh    = 1,
        flag             = 'all',
        gain_key         = args.gain,
        logger           = logger,
        output_dir       = str(outdir),
        phot_apertures   = apertures,
        ref_image        = args.ref_image)

    phot_science = extraction.extract_sources(
        args.fits,
        analysis_thresh  = args.ana_thresh,
        assoc_file       = str(science_xy),
        assoc_cols       = "1,2",
        assoc_radius     = args.host_offset / pix_scale,
        back_size        = args.back_size,
        back_filtersize  = args.back_filtersize,
        deblend_nthresh  = args.deblend_nthresh,
        deblend_mincont  = args.deblend_mincont,
        detect_thresh    = args.det_thresh,
        flag             = 'science',
        gain_key         = args.gain,
        logger           = logger,
        output_dir       = str(outdir),
        phot_apertures   = apertures,
        ref_image        = args.ref_image)

    print(phot_science)

    phot_all     = extraction.postprocess_catalog(phot_all)
    phot_science = extraction.postprocess_catalog(phot_science, verbose=False)

    # Forced photometry if science object not detected
    if not len(phot_science):
        msg = (f'No detection at RA={args.ra}, Dec={args.dec} within '
               f'{args.host_offset:.1f} arcsec. Performing forced photometry.')
        print(bcolors.FAIL + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        with fits.open(args.fits) as hdu:
            hdu_data = hdu[0].data

        forced_phot = extraction.aperture_photometry(
            hdu_data,
            [x_exp, y_exp],
            apertures,
            1.2 * apertures,
            2.0 * apertures,
            1.0,
            gain      = 10,
            zeropoint = np.array(summary_zeropoint['ZP'][2:]))
    else:
        forced_phot = []

    # Apply zeropoint calibration
    for key in [k for k in phot_science.colnames if 'MAG_' in k]:
        zp_match = summary_zeropoint['ZP'][summary_zeropoint['METHOD'] == key]
        if not len(zp_match):
            continue
        zp_val  = float(zp_match[0])
        erp_val = float(summary_zeropoint['ZP_ERRP'][
            summary_zeropoint['METHOD'] == key][0])
        erm_val = float(summary_zeropoint['ZP_ERRM'][
            summary_zeropoint['METHOD'] == key][0])

        errp_key = key.replace('MAG_', 'MAGERRP_')
        errm_key = key.replace('MAG_', 'MAGERRM_')

        if not len(forced_phot):
            phot_science[key]     += zp_val
            phot_science[errp_key] = np.sqrt(phot_science[errp_key]**2 + erp_val**2)
            phot_science[errm_key] = np.sqrt(phot_science[errm_key]**2 + erm_val**2)
            for k in (key, errp_key, errm_key):
                phot_science[k].format = '7.3f'

        phot_all[key]     += zp_val
        phot_all[errp_key] = np.sqrt(phot_all[errp_key]**2 + erp_val**2)
        phot_all[errm_key] = np.sqrt(phot_all[errm_key]**2 + erm_val**2)
        for k in (key, errp_key, errm_key):
            phot_all[k].format = '7.3f'

    # -----------------------------------------------------------------------
    # Step 5: Summaries
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 5: Summaries\n' + bcolors.ENDC)

    j = 0
    for i in range(len(summary_zeropoint['METHOD'])):
        if 'APER' in summary_zeropoint['METHOD'][i]:
            summary_zeropoint['r(FWHM)'][i] = args.ap_diam[j]
            j += 1

    summary_zeropoint['diam(px)']     = 2.0 * summary_zeropoint['r(FWHM)'] * FWHM_median
    summary_zeropoint['diam(arcsec)'] = summary_zeropoint['diam(px)'] * pix_scale
    summary_zeropoint['MAG_3UL_GLOB'] = np.nan * np.ones(len(summary_zeropoint))
    summary_zeropoint['AP_cor']       = (summary_zeropoint['ZP'][-1]
                                         - summary_zeropoint['ZP'])

    for key in summary_zeropoint['METHOD']:
        if key not in phot_all.colnames:
            continue
        errkw = 'MAGERR_' + key.split('_', 1)[1]
        if errkw not in phot_all.colnames:
            continue
        mags    = np.asarray(phot_all[key], dtype=float)
        errs    = np.asarray(phot_all[errkw], dtype=float)
        ul_mask = (errs >= 0.3) & (errs <= 0.34)
        if ul_mask.any():
            summary_zeropoint['MAG_3UL_GLOB'][
                summary_zeropoint['METHOD'] == key] = float(mags[ul_mask].min())

    for k in ('ZP', 'ZP_ERRP', 'ZP_ERRM', 'MAG_3UL_GLOB', 'AP_cor',
               'diam(px)', 'diam(arcsec)'):
        summary_zeropoint[k].format = '.3f'

    print('\n')
    print(bcolors.OKGREEN + '\nZeropoint\n' + bcolors.ENDC)
    logger.info('Zeropoint summary')
    logger.info(str(summary_zeropoint))
    print(summary_zeropoint)

    print(bcolors.OKGREEN + '\nScience photometry\n' + bcolors.ENDC)

    summary_science = calibration.make_scicat(
        args.fits,
        object_properties,
        phot_science,
        forced_phot,
        summary_zeropoint,
        args.host_offset,
        logger)
    logger.info(str(summary_science))
    summary_science.pprint(max_lines=-1)

    # -----------------------------------------------------------------------
    # Step 6: Poststamp
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 6: Poststamp\n' + bcolors.ENDC)

    obs_x_col = summary_science['VALUE'][
        summary_science['PROPERTY'] == 'XWIN_IMAGE_OBS']
    obs_y_col = summary_science['VALUE'][
        summary_science['PROPERTY'] == 'YWIN_IMAGE_OBS']

    if (len(phot_science) > 0
            and 'DISTANCE (arcsec)' in phot_science.colnames
            and float(phot_science['DISTANCE (arcsec)'][0]) <= args.host_offset
            and len(obs_x_col)):
        coord_obs = (float(obs_x_col[0]), float(obs_y_col[0]))
    else:
        coord_obs = (None, None)

    calibration.make_poststamp(
        args.fits, (x_exp, y_exp), coord_obs, output_dir=str(outdir))

    # -----------------------------------------------------------------------
    # Step 7: Save outputs
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 7: Save to file\n' + bcolors.ENDC)

    stem = Path(args.fits).stem
    ascii.write(summary_science,   outdir / (stem + '_phot.log'),         overwrite=True)
    ascii.write(summary_zeropoint, outdir / (stem + '_zp.log'),           overwrite=True)
    ascii.write(phot_all,          outdir / (stem + '_all_abs_cal.phot'), overwrite=True)

    # -----------------------------------------------------------------------
    # Step 8: Cleanup
    # -----------------------------------------------------------------------
    if not args.keeptemp:
        print(bcolors.HEADER + bcolors.BOLD + '\nStep 8: Remove temporary files\n'
              + bcolors.ENDC)
        for pattern in [
            str(outdir / f'check_{Path(args.fits).name}'),
            *glob.glob(str(outdir / (stem + '*refreg*'))),
            *glob.glob(str(outdir / (stem + '*refcat*'))),
            *glob.glob(str(outdir / (stem + '*_loc_*'))),
            *glob.glob(str(outdir / (stem + '_ref_star.*'))),
            str(outdir / (stem + '_xy.cat')),
            str(outdir / (stem + '_all.phot')),
            str(Path('check_' + Path(args.fits).name)),
        ]:
            p = Path(pattern)
            if p.exists():
                try:
                    p.unlink()
                except OSError:
                    pass

    if not args.auto and not args.mag_stdbright and not args.mag_stdfaint:
        plt.show()

    plt.close('all')


if __name__ == '__main__':
    main(sys.argv[1:])
