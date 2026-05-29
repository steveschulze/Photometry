#!/usr/bin/env python
"""
photometry_hst.py — Aperture photometry pipeline for HST drizzled images.

Workflow
--------
1. Verify target in image footprint.
2. Detect sources with sep; optionally centroid on the nearest detection.
3. Estimate background RMS.
4. Aperture photometry at user-supplied diameters.
5. Curve-of-growth analysis.
6. Diagnostic aperture cutout.
7. Write output tables.

Usage
-----
    python photometry_hst.py --ra 11:33:41.550 --dec 00:43:33.50 \\
        --fits SN2015bn_F625W_drc.fits

Run ``python photometry_hst.py --help`` for all options.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import   argparse
import   glob
import   logging
import   sys
from     pathlib import Path

import   numpy as np
from     astropy import table, time
from     astropy.io import ascii, fits
from     matplotlib import pylab as plt

import   extraction
import   fits_tools
import   hst_routines
from     utils import bcolors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    p = argparse.ArgumentParser(
        description='Aperture photometry for HST drizzled images.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--ra',   type=str, required=True,
                   help='RA(J2000) in HMS or decimal degrees.')
    p.add_argument('--dec',  type=str, required=True,
                   help='Dec(J2000) in HMS or decimal degrees.')
    p.add_argument('--fits', type=str, required=True,
                   help='HST drizzled FITS filename.')

    p.add_argument('--ana-thresh',         type=float, default=2,
                   help='Analysis threshold (sigma).')
    p.add_argument('--det-thresh',         type=float, default=2,
                   help='Detection threshold (sigma).')
    p.add_argument('--ap-diam',            type=float, nargs='+',
                   default=[0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00],
                   help='Aperture diameters in arcsec.')
    p.add_argument('--ap-inner-annulus',   type=float, default=1.25,
                   help='Inner annulus scaling factor (× aperture diameter).')
    p.add_argument('--ap-outer-annulus',   type=float, default=2.5,
                   help='Outer annulus scaling factor (× aperture diameter).')

    p.add_argument('--auto',     action='store_true', default=False)
    p.add_argument('--bw',       action='store_true', default=False)
    p.add_argument('--centroid', action='store_true', default=False,
                   help='Centroid on the nearest detected source.')
    p.add_argument('--keeptemp', action='store_true', default=False)
    p.add_argument('--loglevel', type=str, default='INFO')
    p.add_argument('--outdir',   type=str, default='results/')
    p.add_argument('--sex-loglevel', type=str, default='WARNING')
    p.add_argument('--tol',      type=float, default=2,
                   help='Centroid search radius (arcsec).')

    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    """Run the HST aperture photometry pipeline."""
    parser = get_parser()
    args   = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Logger
    log_path = Path(args.fits).with_suffix('.log')
    logger   = logging.getLogger('photometry_hst')
    logger.setLevel(args.loglevel)
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(args.loglevel)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    # -----------------------------------------------------------------------
    # Step 1: Administration
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 1: Administration\n' + bcolors.ENDC)
    logger.info('Output directory: %s', outdir)

    # Convert coordinates
    ra_dd, dec_dd = fits_tools.convert_hms_dd(args.ra, args.dec)

    print(bcolors.OKGREEN + '\nIs the target in the image footprint?' + bcolors.ENDC)
    result = fits_tools.sky2xy(args.fits, ra=ra_dd, dec=dec_dd)
    if result is False:
        msg = (f'Target (RA={args.ra}, Dec={args.dec}) is outside the '
               f'footprint of {args.fits}.')
        print(bcolors.FAIL + msg + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    x_exp, y_exp = result
    x_obs, y_obs = x_exp, y_exp
    print(f'Coordinates = {x_exp:.1f}, {y_exp:.1f}')
    print('Yes.')
    logger.info('Target in footprint: OK')

    science_xy = str(outdir / (Path(args.fits).stem + '_xy.cat'))
    ascii.write(np.array([[x_obs, y_obs]]),
                science_xy, format='no_header', overwrite=True)

    pix2arcsec = fits_tools.pix2arcsec(args.fits)

    # -----------------------------------------------------------------------
    # Step 2: Source catalogue
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 2: Generate source catalogue\n'
          + bcolors.ENDC)
    logger.info('Step 2: source catalogue')

    sources = extraction.extract_sources(
        args.fits,
        analysis_thresh = args.ana_thresh,
        detect_thresh   = args.det_thresh,
        flag            = 'centroid',
        gain_key        = 'CCDGAIN',
        logger          = logger,
        output_dir      = str(outdir),
        phot_apertures  = np.array([1.0]))

    if not len(sources):
        msg = 'No sources detected. Adjust detection thresholds.'
        print(bcolors.FAIL + f'\n{msg}\n' + bcolors.ENDC)
        logger.error(msg)
        sys.exit(1)

    # Distance to expected position
    sources['DISTANCE'] = (
        np.hypot(sources['XWIN_IMAGE'] - x_exp,
                 sources['YWIN_IMAGE'] - y_exp) * pix2arcsec)
    sources['DISTANCE'].unit   = 'arcsec'
    sources['DISTANCE'].format = '.2f'

    # -----------------------------------------------------------------------
    # Optional centroiding
    # -----------------------------------------------------------------------
    if args.centroid:
        nearby = sources[sources['DISTANCE'] <= args.tol]
        if not len(nearby):
            msg = (f'No source within {args.tol:.1f}" of the target. '
                   'Try a larger --tol or lower thresholds.')
            print(bcolors.FAIL + f'\n{msg}\n' + bcolors.ENDC)
            logger.error(msg)
            sys.exit(1)

        nearby.sort('DISTANCE')
        x_obs, y_obs = float(nearby['XWIN_IMAGE'][0]), float(nearby['YWIN_IMAGE'][0])
        ascii.write(np.array([[x_obs, y_obs]]),
                    science_xy, format='no_header', overwrite=True)
        print(bcolors.OKGREEN
              + f'Centroid: ({x_obs:.2f}, {y_obs:.2f})'
              + bcolors.ENDC)
        logger.info('Centroid: (%.2f, %.2f)', x_obs, y_obs)

    # -----------------------------------------------------------------------
    # Step 3: Background RMS
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD
          + '\nStep 3: Background estimation\n' + bcolors.ENDC)
    logger.info('Step 3: background')

    # Load image data for background estimate
    with fits.open(args.fits) as hdulist:
        hdu_header = hdulist[0].header.copy()
        if len(hdulist) > 1:
            try:
                hdulist[1].data.shape
                hdu_header += hdulist[1].header
                hdu_image   = hdulist[1].data
            except (AttributeError, TypeError):
                hdu_header += hdulist[1].header
                hdu_image   = hdulist[0].data
        else:
            hdu_image = hdulist[0].data

    image_rms = extraction.background(hdu_image)

    # -----------------------------------------------------------------------
    # Step 4: Aperture photometry
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 4: Aperture photometry\n'
          + bcolors.ENDC)
    logger.info('Step 4: photometry')

    apertures    = np.array(args.ap_diam)
    inner_annuli = args.ap_inner_annulus * apertures
    outer_annuli = args.ap_outer_annulus * apertures

    if args.ap_inner_annulus < 1:
        print(bcolors.WARNING
              + 'Sky annulus overlaps source aperture. '
              'Check --ap-inner-annulus.'
              + bcolors.ENDC)

    photometry_tbl = hst_routines.hst_aperture_photometry(
        args.fits,
        np.array([x_obs, y_obs]),
        apertures    / 2.0,
        inner_annuli / 2.0,
        outer_annuli / 2.0,
        pix2arcsec,
        image_rms)

    ascii.write(photometry_tbl,
                str(outdir / Path(args.fits).with_suffix('.mag').name),
                overwrite=True)

    # -----------------------------------------------------------------------
    # Step 5: Curve of growth
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 5: Curve of growth\n'
          + bcolors.ENDC)
    logger.info('Step 5: curve of growth')

    hst_routines.hst_cog(
        args.fits,
        np.array([x_obs, y_obs]),
        args.ap_inner_annulus,
        args.ap_outer_annulus,
        pix2arcsec,
        image_rms,
        str(outdir))

    # -----------------------------------------------------------------------
    # Step 6: Diagnostic cutout
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 6: Diagnostic cutout\n'
          + bcolors.ENDC)
    logger.info('Step 6: cutout')

    hst_routines.hst_make_cutout(
        args.fits,
        [x_obs, y_obs],
        [x_exp, y_exp],
        apertures    / 2.0,
        inner_annuli / 2.0,
        outer_annuli / 2.0,
        pix2arcsec,
        str(outdir))

    # -----------------------------------------------------------------------
    # Step 7: Output tables
    # -----------------------------------------------------------------------
    print(bcolors.HEADER + bcolors.BOLD + '\nStep 7: Prepare output\n'
          + bcolors.ENDC)
    logger.info('Step 7: output tables')

    n = len(apertures)

    # Zeropoint table
    zp_tbl              = table.Table()
    zp_tbl['METHOD']    = [f'MAG_APER_{i}' for i in range(n)]
    zp_tbl['ZP']        = [f'{float(photometry_tbl[f"ZP_APER_{i}"][0]):.3f}'
                            for i in range(n)]
    zp_tbl['ZP_ERRP']   = -99.0
    zp_tbl['ZP_ERRM']   = -99.0
    zp_tbl['NUMBER']    = -99
    zp_tbl['r(FWHM)']   = -99.0
    zp_tbl['d(px)']     = apertures / pix2arcsec
    zp_tbl['d(arcsec)'] = apertures
    zp_tbl['MAG_3UL_GLOB'] = [
        float(-2.5 * np.log10(3 * photometry_tbl[f'BKG_NOISE_{i}'][0])
               + photometry_tbl[f'ZP_APER_{i}'][0])
        for i in range(n)]

    for key in ('d(px)', 'd(arcsec)', 'MAG_3UL_GLOB'):
        zp_tbl[key].format = '.3f'

    # Science table
    sci_tbl = table.Table(
        names=('PROPERTY', 'VALUE', 'ERROR+', 'ERROR-', 'COMMENT'),
        dtype=('U100', 'f', 'f', 'f', 'U100'))
    sci_tbl.add_row(['FILENAME', np.nan, np.nan, np.nan, args.fits])

    for key, fallback in (('DATE-OBS', None), ('EXPTIME', None)):
        try:
            if key == 'DATE-OBS':
                t     = time.Time(hdu_header['DATE-OBS'], format='isot',
                                  scale='utc')
                sci_tbl.add_row(['DATE-OBS', np.nan, np.nan, np.nan, t.isot])
                sci_tbl.add_row(['MJD',      np.nan, np.nan, np.nan,
                                  str(round(t.mjd, 7))])
            else:
                sci_tbl.add_row([key, np.nan, np.nan, np.nan,
                                  str(round(hdu_header[key], 2))])
        except Exception:
            sci_tbl.add_row([key, np.nan, np.nan, np.nan, '...'])
            if key == 'DATE-OBS':
                sci_tbl.add_row(['MJD', np.nan, np.nan, np.nan, '...'])

    sci_tbl.add_row(['NCOMBINE', np.nan, np.nan, np.nan, '1'])
    sci_tbl.add_row(['RA',  ra_dd,  np.nan, np.nan, 'degree'])
    sci_tbl.add_row(['DEC', dec_dd, np.nan, np.nan, 'degree'])
    sci_tbl.add_row(['X_IMAGE_EXP', x_exp, np.nan, np.nan, 'px'])
    sci_tbl.add_row(['Y_IMAGE_EXP', y_exp, np.nan, np.nan, 'px'])
    sci_tbl.add_row(['X_IMAGE_OBS', round(x_obs, 3), np.nan, np.nan, 'px'])
    sci_tbl.add_row(['Y_IMAGE_OBS', round(y_obs, 3), np.nan, np.nan, 'px'])

    dist_arcsec = float(np.hypot(x_obs - x_exp, y_obs - y_exp) * pix2arcsec)
    sci_tbl.add_row(['DISTANCE (px)',
                      round(np.hypot(x_obs - x_exp, y_obs - y_exp), 3),
                      np.nan, np.nan, 'px'])
    sci_tbl.add_row(['DISTANCE (arcsec)', round(dist_arcsec, 3),
                      np.nan, np.nan, 'arcsec'])

    for i in range(n):
        fnu     = float(photometry_tbl[f'FNU_APER_{i}'][0])
        fnu_err = float(photometry_tbl[f'FNUERR_APER_{i}'][0])
        zp_i    = float(photometry_tbl[f'ZP_APER_{i}'][0])

        sci_tbl.add_row([f'FNU_APER_{i}',
                          f'{fnu:.3e}', f'{fnu_err:.3e}', f'{fnu_err:.3e}',
                          'microJy'])
        if fnu > 0:
            sci_tbl.add_row([f'MAG_APER_{i}',
                              round(float(photometry_tbl[f'MAG_APER_{i}'][0]),  3),
                              round(float(photometry_tbl[f'MAGERRP_APER_{i}'][0]), 3),
                              round(float(photometry_tbl[f'MAGERRM_APER_{i}'][0]), 3),
                              'mag'])
        else:
            sci_tbl.add_row([f'MAG_APER_{i}',
                              round(-2.5 * np.log10(3 * fnu_err) + 23.9, 3),
                              np.nan, np.nan, 'mag'])

        sci_tbl.add_row([f'MAG_APER_{i}_2sigma',
                          round(-2.5 * np.log10(2 * fnu_err) + 23.9, 3),
                          np.nan, np.nan, 'mag'])
        sci_tbl.add_row([f'MAG_APER_{i}_3sigma',
                          round(-2.5 * np.log10(3 * fnu_err) + 23.9, 3),
                          np.nan, np.nan, 'mag'])

    # Print and save
    print(bcolors.OKGREEN + '\nZeropoint\n' + bcolors.ENDC)
    print(zp_tbl)
    print(bcolors.OKGREEN + '\nScience\n' + bcolors.ENDC)
    sci_tbl.pprint(max_lines=-1)

    stem = Path(args.fits).stem
    ascii.write(zp_tbl,  outdir / (stem + '_zp.log'),   overwrite=True)
    ascii.write(sci_tbl, outdir / (stem + '_phot.log'), overwrite=True)

    # -----------------------------------------------------------------------
    # Step 8: Cleanup
    # -----------------------------------------------------------------------
    if not args.keeptemp:
        print(bcolors.HEADER + bcolors.BOLD + '\nStep 8: Cleanup\n' + bcolors.ENDC)
        for p in [
            Path('check_' + Path(args.fits).name),
            outdir / (stem + '_xy.cat'),
            outdir / (stem + '_centroid.log'),
        ]:
            if p.exists():
                try:
                    p.unlink()
                except OSError:
                    pass

    if not args.auto:
        plt.show()
    plt.close('all')


if __name__ == '__main__':
    main(sys.argv[1:])
