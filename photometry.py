#!/usr/bin/env python
"""
photometry.py — Aperture photometry pipeline for ground-based imaging.

Workflow
--------
1. Verify that the target lies in the image footprint (WCS check).
2. Retrieve or read a reference photometric catalogue.
3. Build a local comparison-star sequence (SExtractor + cross-match).
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
from     astropy import table, time
from     astropy.io import ascii, fits
from     astropy import units as u
import   cat_tools
import   copy
import   fits_tools
import   logging
from     matplotlib import pylab as plt
import   numpy as np
import   numpy.lib.recfunctions as rfn
import   os
from     pathlib import Path
import   phot_routines
from     plotsettings import *
import   sys

# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def get_parser():
    """Build the argument parser for the photometry pipeline.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(
        description='Aperture photometry pipeline for ground-based images.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Target coordinates
    p.add_argument('--ra',  type=str, required=True,
                   help='RA(J2000) in HMS or decimal degrees. '
                        'For negative Dec write \\" -12:20:20.2\\"')
    p.add_argument('--dec', type=str, required=True,
                   help='Dec(J2000) in DMS or decimal degrees.')
    p.add_argument('--fits', type=str, required=True,
                   help='Science FITS filename.')
    p.add_argument('--host-offset', type=float, default=5,
                   help='Detection search radius around target (arcsec).')

    # Reference catalogue
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

    # Photometry
    p.add_argument('--ana-thresh',       type=float, default=1,
                   help='Analysis threshold (sigma).')
    p.add_argument('--ap-diam',          type=float, default=[1., 1.5, 2., 3.],
                   nargs='+',
                   help='Aperture diameters as multiples of FWHM.')
    p.add_argument('--ap-diam-ul',       type=float, default=2, nargs='+',
                   help='Aperture diameter for forced-phot upper limit.')
    p.add_argument('--det-thresh',       type=float, default=1,
                   help='Detection threshold (sigma).')
    p.add_argument('--gain',             type=str,   default=None,
                   help='Header keyword for detector gain (e-/ADU).')
    p.add_argument('--back-size',        type=int,   default=64,
                   help='Background mesh size (pixels).')
    p.add_argument('--back-filtersize',  type=int,   default=3,
                   help='Background filter size.')
    p.add_argument('--deblend-nthresh',  type=int,   default=64,
                   help='Deblending sub-thresholds.')
    p.add_argument('--deblend-mincont',  type=float, default=0.0001,
                   help='Minimum deblending contrast.')
    p.add_argument('--noflags',          action='store_true', default=False,
                   help='Ignore sep FLAGS when filtering stars.')

    # Zeropoint
    p.add_argument('--mag-cut',       type=float, default=12.,
                   help='Exclude stars brighter than this magnitude.')
    p.add_argument('--mag-stdfaint',  type=str,   default='0',
                   help='Faint magnitude limit for comparison stars.')
    p.add_argument('--mag-stdbright', type=str,   default='0',
                   help='Bright magnitude limit for comparison stars.')
    p.add_argument('--maxstars',      type=int,   default=200,
                   help='Maximum comparison stars.')

    # Control
    p.add_argument('--auto',      action='store_true', default=False,
                   help='Non-interactive mode (auto magnitude cuts).')
    p.add_argument('--bw',        action='store_true', default=False,
                   help='Monochrome terminal output.')
    p.add_argument('--keeptemp',  action='store_true', default=False,
                   help='Keep temporary files after completion.')
    p.add_argument('--loglevel',  type=str, default='INFO',
                   help='Log level (DEBUG/INFO/WARNING/ERROR/CRITICAL).')
    p.add_argument('--outdir',    type=str, default='results/',
                   help='Output directory.')
    p.add_argument('--sex-loglevel', type=str, default='WARNING',
                   help='Logging level passed to sep (currently unused).')
    p.add_argument('--tol',       type=float, default=1,
                   help='Cross-matching tolerance (arcsec).')

    return p


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(argv=None):
    """Run the aperture photometry pipeline.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments.  ``sys.argv[1:]`` is used when None.
    """
    parser = get_parser()
    args   = parser.parse_args(argv)

    args.mag_stdfaint  = float(args.mag_stdfaint)
    args.mag_stdbright = float(args.mag_stdbright)

    # Terminal colour output
    if args.bw:
        class bcolors:
            HEADER = OKBLUE = OKGREEN = WARNING = FAIL = ENDC = BOLD = UNDERLINE = ''
    else:
        from misc import bcolors

    # Logger setup
    logger = logging.getLogger('photometry')
    logger.setLevel(args.loglevel)
    log_path = Path(args.fits).with_suffix('.log')
    fh       = logging.FileHandler(log_path, mode='w')
    fh.setLevel(args.loglevel)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    # -----------------------------------------------------------------------
    # Step 1: Administration
    # -----------------------------------------------------------------------
    msg = 'Step 1: Administration'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    outdir = Path(args.outdir)
    if not str(outdir).endswith('/'):
        outdir = Path(str(outdir) + '/')
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info(f'Output directory: {outdir}')

    # Reconstruct command string for logging
    cmd = 'photometry.py '
    for key in sorted(vars(args)):
        val = vars(args)[key]
        if key in ('auto', 'bw', 'noflags') and val:
            cmd += f'--{key} '
        elif key == 'ref_image' and val:
            cmd += f'--{key.replace("_","-")} {val} '
        elif key in ('ref_cat', 'ref_filter', 'ref_file'):
            if vars(args)['ref_file'] is not None:
                if '--ref-file' not in cmd:
                    cmd += f'--ref-file {vars(args)["ref_file"]} '
            elif key != 'ref_file' and val is not None:
                cmd += f'--{key.replace("_","-")} {val} '
        elif key == 'ap_diam':
            cmd += f'--ap-diam {" ".join(str(x) for x in val)} '
        elif key == 'dec':
            v = f'"{val}"' if ('-' in str(val) or '+' in str(val)) else val
            cmd += f'--dec {v} '
        else:
            if val not in (None, '', False):
                cmd += f'--{key.replace("_","-")} {val} '

    print(bcolors.OKGREEN + '\nCommand:\n' + cmd + bcolors.ENDC)
    logger.info(f'Command: {cmd}')

    ra_dd, dec_dd = fits_tools.convert_hms_dd(args.ra, args.dec)

    object_properties = table.Table()
    object_properties['OBJECT'] = [Path(args.fits).stem]
    object_properties['RA_HMS'] = [args.ra]
    object_properties['DEC_DMS'] = [args.dec]
    object_properties['RA']     = [ra_dd]
    object_properties['DEC']    = [dec_dd]

    # Verify target is in image footprint
    print(bcolors.OKGREEN + '\nIs the target in the image footprint?' + bcolors.ENDC)
    result = fits_tools.sky2xy(args.fits, RA=ra_dd, DEC=dec_dd)
    if result is False:
        msg = (f'Target (RA={args.ra}, Dec={args.dec}) is outside the image '
               f'footprint of {args.fits}.')
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    x_exp, y_exp = result
    object_properties['X_EXP'] = [x_exp]
    object_properties['Y_EXP'] = [y_exp]

    filename_science_xy = outdir / (object_properties['OBJECT'][0] + '_xy.cat')
    ascii.write(np.array([[x_exp, y_exp]]),
                str(filename_science_xy), format='no_header', overwrite=True)
    print('Yes.')
    logger.info('Target in footprint: OK')

    # -----------------------------------------------------------------------
    # Step 2: Reference catalogue
    # -----------------------------------------------------------------------
    msg = 'Step 2: Flux calibration'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    if args.ref_file is not None:
        msg = f'Using catalogue file: {args.ref_file}'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        ref_cat = ascii.read(args.ref_file,
                             names=('RA', 'DEC', 'MAG', 'MAGERR'))
        object_properties['PHOTCAL'] = [args.ref_file]

        base = (Path(args.ref_file).stem)
        filename_stars = outdir / (object_properties['OBJECT'][0]
                                   + '_' + base + '.cat')
        ascii.write(ref_cat, str(filename_stars).replace('.cat', '_refcat.cat'),
                    overwrite=True, format='no_header')
        logger.info(f'Reference catalogue copy: {filename_stars}')

    else:
        msg = 'Retrieving catalogue from Vizier'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        object_properties['PHOTCAL'] = [args.ref_cat + '/' + args.ref_filter]
        object_properties['CAT']     = [args.ref_cat]
        object_properties['FILTER']  = [args.ref_filter]

        filename_stars = outdir / (object_properties['OBJECT'][0]
                                   + '_' + args.ref_cat
                                   + '_' + args.ref_filter + '.cat')

        cat_tools.retrieve_photcat(
            object_properties, args.ref_cat, phot_routines.catalog_prop,
            FILENAME=str(filename_stars),
            RADIUS=args.ref_radius * u.arcmin,
            OUTDIR=str(outdir))

        refcat_path = str(filename_stars).replace('.cat', '_refcat.cat')
        import shutil
        shutil.copy(str(filename_stars), refcat_path)
        logger.info(f'Reference catalogue copy: {refcat_path}')

    # Build local sequence
    msg = 'Building the local sequence'
    print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
    logger.info(msg)

    refcat_path = str(filename_stars).replace('.cat', '_refcat.cat')
    ref_cat     = ascii.read(refcat_path, names=('RA', 'DEC', 'MAG', 'MAGERR'))

    # Apply magnitude cut
    ref_cat = ref_cat[ref_cat['MAG'] > args.mag_cut]
    logger.info(f'Removing stars brighter than {args.mag_cut:.2f} mag')
    ascii.write(ref_cat, refcat_path, overwrite=True, format='no_header')

    # Remove stars outside image footprint (vectorised)
    refcat_xy_path = refcat_path.replace('_refcat.cat', '_refcat_xy.cat')
    catalog_clean  = fits_tools.sky2xy(args.fits, CAT=refcat_path)

    if len(catalog_clean) == 0:
        msg = 'No reference stars in the image footprint.'
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    ascii.write(catalog_clean, refcat_xy_path, overwrite=True, format='no_header')
    print(bcolors.OKGREEN + f'\n{len(catalog_clean)} reference stars in footprint.'
          + bcolors.ENDC)

    # -----------------------------------------------------------------------
    # Run sep on reference stars (local sequence detection)
    # -----------------------------------------------------------------------
    print(bcolors.OKGREEN + '\nSelecting stars for local sequence' + bcolors.ENDC)

    pix_scale = fits_tools.pix2arcsec(args.fits)

    ref_stars = phot_routines.sextractor_photometry(
        ANALYSIS_THRESH  = 3,
        ASSOC_NAME       = refcat_xy_path,
        ASSOC_PARAMS     = "1,2",
        ASSOC_RADIUS     = args.tol / pix_scale,
        DETECT_THRESH    = 3,
        FITS             = args.fits,
        FLAG             = 'loc_seq',
        GAIN             = args.gain,
        LOGGER           = logger,
        LOGLEVEL         = args.sex_loglevel,
        PATH             = str(outdir),
        PHOT_APERTURES   = np.array([10.]),
        REF_FILE         = args.ref_image)

    if len(ref_stars) == 0:
        msg = 'No reference stars detected. Check --tol, image WCS, and catalog.'
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    # Quality cuts
    ref_stars = ref_stars[
        (ref_stars['MAGERR_APER'] > 0.0002)
        & (ref_stars['MAGERR_AUTO'] > 0.0002)]
    logger.info(f'{len(ref_stars)} stars identified for local sequence')

    if args.mag_stdbright == 0 and args.mag_stdfaint == 0:
        logger.info('Pruning: flags=0, FWHM IQR clipping, max stars limit')

        if not args.noflags:
            ref_stars = ref_stars[ref_stars['FLAGS'] == 0]

        # PSF size clipping (IQR fence)
        fwhm_col = 'FWHM_IMAGE'
        if np.all(np.array(ref_stars[fwhm_col]) == 0):
            print(bcolors.WARNING
                  + 'All FWHM_IMAGE = 0; using FLUX_RADIUS instead.'
                  + bcolors.ENDC)
            fwhm_vals = np.array(ref_stars['FLUX_RADIUS']) * 2. / 1.1
        else:
            fwhm_vals = np.array(ref_stars[fwhm_col])

        p25, p75 = np.percentile(fwhm_vals[fwhm_vals > 0], [25, 75])
        iqr      = p75 - p25
        fwhm_ok  = fwhm_vals < (p75 + 1.5 * iqr)
        ref_stars = ref_stars[fwhm_ok]

        if len(ref_stars) > args.maxstars:
            print(bcolors.WARNING
                  + f'Truncating to brightest {args.maxstars} stars.'
                  + bcolors.ENDC)
            ref_stars.sort('MAGERR_AUTO')
            ref_stars = ref_stars[:args.maxstars]

        logger.info(f'{len(ref_stars)} stars remain after pruning')

        # Cross-match sep catalog with reference catalog
        msg = 'Cross-matching catalogues'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        ref_stars_keys = [
            'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE',
            'MAG_AUTO', 'MAGERR_AUTO', 'MAG_PETRO', 'MAGERR_PETRO',
            'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_PETRO', 'FLUXERR_PETRO',
            'FWHM_IMAGE', 'FWHM_WORLD', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
            'KRON_RADIUS', 'FLUX_RADIUS', 'FLAGS', 'VECTOR_ASSOC',
            'NUMBER_ASSOC', 'MAG_APER', 'MAGERR_APER',
            'FLUX_APER', 'FLUXERR_APER']
        ref_stars_keys = [k for k in ref_stars_keys
                          if k in ref_stars.colnames]
        ref_cat_keys   = list(ref_cat.colnames)

        # numpy 2.x-safe conversion (replaces .view((float, n)))
        ref_stars_arr = rfn.structured_to_unstructured(
            np.asarray(ref_stars[ref_stars_keys]), dtype=float)
        ref_cat_arr   = rfn.structured_to_unstructured(
            np.asarray(ref_cat[ref_cat_keys]), dtype=float)

        matched_arr = cat_tools.wrapper_crossmatch(
            ref_stars_arr, ref_cat_arr, args.tol)

        if len(matched_arr) == 0:
            msg = (f'No cross-matches found within {args.tol} arcsec. '
                   'Check --tol and WCS accuracy.')
            print(bcolors.FAIL + msg + bcolors.ENDC)
            logger.error(msg)
            sys.exit(1)

        matched_standard = table.Table(
            matched_arr,
            names=ref_cat_keys + ref_stars_keys + ['DIST'])
        matched_standard['DIST'] = matched_standard['DIST'] * 3600.
        matched_standard['DIST'].format = '.3f'

        matched_standard.rename_column('MAG_APER',    'MAG_INS')
        matched_standard.rename_column('MAGERR_APER', 'MAGERR_INS')
        matched_standard.rename_column('MAG',         'MAG_CAT')
        matched_standard.rename_column('MAGERR',      'MAGERR_CAT')
        matched_standard['MAG_ZP_TEMP'] = (matched_standard['MAG_CAT']
                                            - matched_standard['MAG_INS'])
        matched_standard['MAG_ZP_TEMP'].format = '.3f'

        temp = copy.deepcopy(matched_standard)
        temp.sort('MAG_CAT')
        print(temp[['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE',
                     'YWIN_IMAGE', 'FWHM_IMAGE', 'MAG_CAT', 'MAGERR_CAT',
                     'MAG_INS', 'MAGERR_INS', 'MAG_ZP_TEMP', 'DIST']])

        local_seq = phot_routines.local_sequence(
            matched_standard,
            AUTO     = args.auto,
            FILENAME = str(filename_stars).replace('.cat', '_loc_cleaned.cat'),
            FITS     = args.fits,
            LOGGER   = logger,
            LOWER    = 5,
            PATH     = str(outdir),
            UPPER    = 90)

        matched_standard = local_seq['CAT']
        logger.info(f'Local sequence stars: {local_seq["NUMSTARS"]}')
        print(bcolors.OKGREEN
              + f'\nFinal local sequence: {local_seq["NUMSTARS"]} stars'
              + bcolors.ENDC)

    else:
        # Manual magnitude range
        logger.info(f'Manual magnitude range: {args.mag_stdbright:.2f} '
                    f'– {args.mag_stdfaint:.2f}')
        ref_stars = ref_stars[
            (ref_stars['MAG_APER'] >= args.mag_stdbright)
            & (ref_stars['MAG_APER'] <= args.mag_stdfaint)]

        if not args.noflags:
            ref_stars = ref_stars[ref_stars['FLAGS'] == 0]

        fwhm_vals = np.array(ref_stars['FWHM_IMAGE'])
        if not np.all(fwhm_vals == 0):
            p25, p75 = np.percentile(fwhm_vals[fwhm_vals > 0], [25, 75])
            iqr      = p75 - p25
            ref_stars = ref_stars[fwhm_vals < p75 + 1.5 * iqr]

        if len(ref_stars) > args.maxstars:
            ref_stars.sort('MAGERR_AUTO')
            ref_stars = ref_stars[:args.maxstars]

        msg = 'Cross-matching catalogues'
        print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        ref_stars_keys = [
            'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE',
            'MAG_AUTO', 'MAGERR_AUTO', 'MAG_PETRO', 'MAGERR_PETRO',
            'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_PETRO', 'FLUXERR_PETRO',
            'FWHM_IMAGE', 'FWHM_WORLD', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
            'KRON_RADIUS', 'FLUX_RADIUS', 'FLAGS', 'VECTOR_ASSOC',
            'NUMBER_ASSOC', 'MAG_APER', 'MAGERR_APER',
            'FLUX_APER', 'FLUXERR_APER']
        ref_stars_keys = [k for k in ref_stars_keys
                          if k in ref_stars.colnames]
        ref_cat_keys   = list(ref_cat.colnames)

        ref_stars_arr = rfn.structured_to_unstructured(
            np.asarray(ref_stars[ref_stars_keys]), dtype=float)
        ref_cat_arr   = rfn.structured_to_unstructured(
            np.asarray(ref_cat[ref_cat_keys]), dtype=float)

        matched_arr = cat_tools.wrapper_crossmatch(
            ref_stars_arr, ref_cat_arr, args.tol)

        if len(matched_arr) == 0:
            msg = f'No cross-matches found within {args.tol} arcsec.'
            print(bcolors.FAIL + msg + bcolors.ENDC)
            logger.error(msg)
            sys.exit(1)

        matched_standard = table.Table(
            matched_arr,
            names=ref_cat_keys + ref_stars_keys + ['DIST'])
        matched_standard['DIST'] = matched_standard['DIST'] * 3600.
        matched_standard['DIST'].format = '.3f'

        matched_standard.rename_column('MAG_APER',    'MAG_INS')
        matched_standard.rename_column('MAGERR_APER', 'MAGERR_INS')
        matched_standard.rename_column('MAG',         'MAG_CAT')
        matched_standard.rename_column('MAGERR',      'MAGERR_CAT')

        local_seq = phot_routines.local_sequence(
            matched_standard,
            AUTO     = args.auto,
            FILENAME = str(filename_stars).replace('.cat', '_loc_cleaned.cat'),
            FITS     = args.fits,
            LOGGER   = logger,
            LOWER    = 0,
            PATH     = str(outdir),
            UPPER    = 100)

        matched_standard = local_seq['CAT']
        logger.info(f'Local sequence stars: {local_seq["NUMSTARS"]}')
        print(bcolors.OKGREEN
              + f'\nFinal local sequence: {local_seq["NUMSTARS"]} stars'
              + bcolors.ENDC)

    # Determine FWHM for aperture scaling
    if args.ref_image == '':
        fwhm_col = np.array(matched_standard['FWHM_IMAGE'])
        if np.all(fwhm_col == 0):
            print(bcolors.WARNING
                  + 'All FWHM_IMAGE = 0; using FLUX_RADIUS.'
                  + bcolors.ENDC)
            FWHM_median = float(np.median(matched_standard['FLUX_RADIUS'])) * 2. / 1.1
        else:
            FWHM_median = float(np.median(fwhm_col[fwhm_col > 0]))
    else:
        FWHM_median = float(p75)  # noqa: use FWHM from PSF clipping

    apertures = 2. * FWHM_median * np.array(args.ap_diam)
    logger.info(f'FWHM = {FWHM_median:.2f} px; apertures = {apertures} px')

    # -----------------------------------------------------------------------
    # Step 3: Zeropoint
    # -----------------------------------------------------------------------
    msg = 'Step 3: Zeropoint calculation'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    print(bcolors.OKGREEN + bcolors.BOLD + '\nRun sep on comparison stars'
          + bcolors.ENDC)
    logger.info('Run sep on comparison stars')

    phot_stars = phot_routines.sextractor_photometry(
        ANALYSIS_THRESH  = 3,
        ASSOC_NAME       = str(filename_stars).replace('.cat', '_loc_cleaned.cat'),
        ASSOC_PARAMS     = "1,2",
        ASSOC_RADIUS     = args.tol / pix_scale,
        DETECT_THRESH    = 3,
        FITS             = args.fits,
        FLAG             = 'ref_star',
        GAIN             = args.gain,
        LOGGER           = logger,
        LOGLEVEL         = args.sex_loglevel,
        PATH             = str(outdir),
        PHOT_APERTURES   = apertures,
        REF_FILE         = args.ref_image)

    print(bcolors.OKGREEN + bcolors.BOLD + '\nCompute zeropoint' + bcolors.ENDC)
    logger.info('Compute zeropoint')

    summary_zeropoint = phot_routines.zeropoint(
        matched_standard, phot_stars,
        NITER     = 1000,
        FITS      = args.fits,
        LOGGER    = logger,
        PATH      = str(outdir),
        TOLERANCE = args.tol)

    summary_zeropoint['r(FWHM)'] = np.nan

    # -----------------------------------------------------------------------
    # Step 4: Science photometry
    # -----------------------------------------------------------------------
    msg = 'Step 4: Aperture photometry'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    phot_all = phot_routines.sextractor_photometry(
        ANALYSIS_THRESH  = 1,
        BACK_SIZE        = args.back_size,
        BACK_FILTERSIZE  = args.back_filtersize,
        DEBLEND_NTHRESH  = args.deblend_nthresh,
        DEBLEND_MINCONT  = args.deblend_mincont,
        DETECT_THRESH    = 1,
        FLAG             = 'all',
        FITS             = args.fits,
        GAIN             = args.gain,
        LOGGER           = logger,
        LOGLEVEL         = args.sex_loglevel,
        PATH             = str(outdir),
        PHOT_APERTURES   = apertures,
        REF_FILE         = args.ref_image)

    phot_science = phot_routines.sextractor_photometry(
        ANALYSIS_THRESH  = args.ana_thresh,
        ASSOC_NAME       = str(filename_science_xy),
        ASSOC_PARAMS     = "1,2",
        ASSOC_RADIUS     = args.host_offset / pix_scale,
        BACK_SIZE        = args.back_size,
        BACK_FILTERSIZE  = args.back_filtersize,
        DEBLEND_NTHRESH  = args.deblend_nthresh,
        DEBLEND_MINCONT  = args.deblend_mincont,
        DETECT_THRESH    = args.det_thresh,
        FITS             = args.fits,
        FLAG             = 'science',
        GAIN             = args.gain,
        LOGGER           = logger,
        LOGLEVEL         = args.sex_loglevel,
        PATH             = str(outdir),
        PHOT_APERTURES   = apertures,
        REF_FILE         = args.ref_image)

    print(phot_science)

    print(bcolors.OKGREEN + bcolors.BOLD
          + '\nPost-process sep output' + bcolors.ENDC)
    logger.info('Post-process sep output')

    phot_all     = phot_routines.sextractor_postprocess(phot_all)
    phot_science = phot_routines.sextractor_postprocess(phot_science,
                                                         PRINTHELP=False)

    # Forced photometry if science object not detected
    if len(phot_science) == 0:
        msg = (f'No detection at RA={args.ra}, Dec={args.dec} within '
               f'{args.host_offset:.1f} arcsec. Performing forced photometry.')
        print(bcolors.FAIL + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
        logger.info(msg)

        with fits.open(args.fits) as hdu:
            hdu_data = hdu[0].data

        forced_phot = phot_routines.aperture_photometry(
            hdu_data,
            [x_exp, y_exp],
            apertures,
            1.2 * apertures,
            2.0 * apertures,
            1,
            GAIN      = 10,
            ZEROPOINT = np.array(summary_zeropoint['ZP'][2:]),
            FA        = 1)
    else:
        forced_phot = []

    # Apply zeropoint calibration
    msg = 'Convert instrumental to calibrated magnitudes'
    print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
    logger.info(msg)

    mag_keys = [k for k in phot_science.colnames if 'MAG_' in k]
    for key in mag_keys:
        zp_match = summary_zeropoint['ZP'][summary_zeropoint['METHOD'] == key]
        if len(zp_match) == 0:
            continue
        zp_val  = float(zp_match[0])
        erp_val = float(summary_zeropoint['ZP_ERRP'][
            summary_zeropoint['METHOD'] == key][0])
        erm_val = float(summary_zeropoint['ZP_ERRM'][
            summary_zeropoint['METHOD'] == key][0])

        errp_key = key.replace('MAG_', 'MAGERRP_')
        errm_key = key.replace('MAG_', 'MAGERRM_')

        if len(forced_phot) == 0:
            phot_science[key]     += zp_val
            phot_science[errp_key] = np.sqrt(
                phot_science[errp_key]**2 + erp_val**2)
            phot_science[errm_key] = np.sqrt(
                phot_science[errm_key]**2 + erm_val**2)
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
    msg = 'Step 5: Summaries'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    # Fill in aperture-size columns
    j = 0
    for i in range(len(summary_zeropoint['METHOD'])):
        if 'APER' in summary_zeropoint['METHOD'][i]:
            summary_zeropoint['r(FWHM)'][i] = args.ap_diam[j]
            j += 1

    summary_zeropoint['diam(px)']    = (2. * summary_zeropoint['r(FWHM)']
                                         * FWHM_median)
    summary_zeropoint['diam(arcsec)'] = (summary_zeropoint['diam(px)']
                                          * pix_scale)

    # 3σ limiting magnitude from global blank-sky distribution
    summary_zeropoint['MAG_3UL_GLOB'] = np.nan * np.ones(len(summary_zeropoint))
    summary_zeropoint['AP_cor'] = (summary_zeropoint['ZP'][-1]
                                    - summary_zeropoint['ZP'])

    for key in summary_zeropoint['METHOD']:
        if key not in phot_all.colnames:
            continue
        errkw = 'MAGERR_' + key.split('_', 1)[1]
        if errkw not in phot_all.colnames:
            continue
        mags     = np.array(phot_all[key], dtype=float)
        errs     = np.array(phot_all[errkw], dtype=float)
        ul_mask  = (errs >= 0.3) & (errs <= 0.34)
        if ul_mask.any():
            summary_zeropoint['MAG_3UL_GLOB'][
                summary_zeropoint['METHOD'] == key] = float(np.min(mags[ul_mask]))

    for k in ('ZP', 'ZP_ERRP', 'ZP_ERRM', 'MAG_3UL_GLOB', 'AP_cor',
               'diam(px)', 'diam(arcsec)'):
        summary_zeropoint[k].format = '.3f'

    print('\n')
    print(bcolors.OKGREEN + '\nZeropoint summary\n' + bcolors.ENDC)
    logger.info('Zeropoint summary')
    logger.info(str(summary_zeropoint))
    print(summary_zeropoint)

    print(bcolors.OKGREEN + '\nScience photometry summary\n' + bcolors.ENDC)

    summary_science = phot_routines.make_scicat(
        args.fits, object_properties,
        phot_science, forced_phot,
        summary_zeropoint, args.host_offset, logger)
    logger.info('Science photometry summary')
    logger.info(str(summary_science))
    summary_science.pprint(max_lines=-1)

    # -----------------------------------------------------------------------
    # Step 6: Poststamp
    # -----------------------------------------------------------------------
    msg = 'Step 6: Poststamp'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    obs_row = summary_science[summary_science['PROPERTY'] == 'XWIN_IMAGE_OBS']
    if (len(phot_science) > 0
            and 'DISTANCE (arcsec)' in phot_science.colnames
            and float(phot_science['DISTANCE (arcsec)'][0]) <= args.host_offset):
        coord_obs = [float(obs_row['VALUE'][0]),
                     float(summary_science['VALUE'][
                         summary_science['PROPERTY'] == 'YWIN_IMAGE_OBS'][0])]
    else:
        coord_obs = [None, None]

    phot_routines.make_poststamp(
        args.fits, [x_exp, y_exp], coord_obs, PATH=str(outdir))

    # -----------------------------------------------------------------------
    # Step 7: Save outputs
    # -----------------------------------------------------------------------
    msg = 'Step 7: Save to file'
    print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
    logger.info(msg)

    fits_stem = Path(args.fits).stem
    ascii.write(summary_science,   outdir / (fits_stem + '_phot.log'),
                overwrite=True)
    ascii.write(summary_zeropoint, outdir / (fits_stem + '_zp.log'),
                overwrite=True)
    ascii.write(phot_all,          outdir / (fits_stem + '_all_abs_cal.phot'),
                overwrite=True)

    # -----------------------------------------------------------------------
    # Step 8: Cleanup
    # -----------------------------------------------------------------------
    if not args.keeptemp:
        msg = 'Step 8: Remove temporary files'
        print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
        logger.info(msg)

        stem = Path(args.fits).stem
        for pattern in (
            'check_' + Path(args.fits).name,
        ):
            p = Path(pattern)
            if p.exists():
                p.unlink()

        for pattern in (
            str(outdir / (stem + '*refreg*')),
            str(outdir / (stem + '*refcat*')),
            str(outdir / (stem + '*_loc_*')),
            str(outdir / (stem + '_ref_star.*')),
            str(outdir / (stem + '_xy.cat')),
            str(outdir / (stem + '_all.phot')),
        ):
            import glob
            for f in glob.glob(pattern):
                try:
                    Path(f).unlink()
                except OSError:
                    pass
    else:
        msg = 'Step 8: Keeping temporary files'
        print(bcolors.HEADER + bcolors.BOLD + f'\n{msg}\n' + bcolors.ENDC)
        logger.info(msg)

    # Show interactive plots if not in auto mode
    if not args.auto and not args.mag_stdbright and not args.mag_stdfaint:
        plt.show()

    plt.close('all')


if __name__ == '__main__':
    main(sys.argv[1:])
