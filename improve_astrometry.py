#!/usr/bin/env python
"""
improve_astrometry.py — Improve the WCS calibration of a FITS image.

Uses astrometry.net's ``solve-field`` command for blind astrometric
calibration using astrometry.net's own built-in ``image2xy`` source
detector (no SExtractor binary required), then converts the SIP distortion
keywords to PV format for compatibility with sep.

Usage
-----
    python improve_astrometry.py --ra 11:33:41.550 --dec 00:43:33.50 \\
        --fits SN2015bn_SDSS_r_wcs.fits
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import    argparse
from      astropy.io import fits
import    fits_tools
from      utils import bcolors
from      pathlib import Path
import    sip_to_pv
import    subprocess
import    sys


def get_parser():
    """Build the argument parser for improve_astrometry.

    Returns
    -------
    argparse.ArgumentParser
    """
    p = argparse.ArgumentParser(
        description=(
            'Improve WCS calibration with astrometry.net '
            'and convert SIP distortion terms to PV format.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--fits',         type=str, required=True,
                   help='Input FITS filename.')
    p.add_argument('--ra',           type=str, required=True,
                   help='RA(J2000) in HMS or decimal degrees.')
    p.add_argument('--dec',          type=str, required=True,
                   help='Dec(J2000) in HMS or decimal degrees.')
    p.add_argument('--radius',       type=float, default=0.125,
                   help='Search radius around field centre (deg).')
    p.add_argument('--downsample',   type=int,   default=2,
                   help='Downsample factor before source extraction.')
    p.add_argument('--small-field',  action='store_true', default=False,
                   help='Set for fields smaller than 5 arcmin.')
    p.add_argument('--tweak-order',  type=int,   default=2,
                   help='Polynomial order for SIP WCS corrections (0 = no SIP).')
    return p


def main(args=None):
    """Run the astrometry improvement pipeline.

    Parameters
    ----------
    args : list of str, optional
        Command-line arguments. ``sys.argv[1:]`` when None.
    """
    parser = get_parser()
    args   = parser.parse_args(args)

    fits_path = Path(args.fits)
    outfile   = fits_path.with_name(fits_path.stem + '_wcs' + fits_path.suffix)

    log_astro   = open('astro.log',      'w')
    log_sip2pv  = open('astro_sip2pv.log', 'w')

    ra_dd, dec_dd = fits_tools.convert_hms_dd(args.ra, args.dec)

    print(bcolors.HEADER + 'Process image: ' + str(fits_path) + bcolors.ENDC)

    # Sanitise FITS headers
    with fits.open(str(fits_path)) as hdulist:
        data   = hdulist[0].data
        header = hdulist[0].header
    fits.writeto(str(fits_path), data, header,
                 overwrite=True, output_verify='silentfix')

    # Build solve-field command.
    # We do NOT pass --source-extractor-config: astrometry.net's built-in
    # image2xy handles source detection without requiring a SExtractor binary.
    # Using source-extractor externally caused segfaults because the relative
    # paths in default.sex (default.conv, default.nnw) break when solve-field
    # spawns source-extractor from a temp directory.
    cmd_base = [
        'solve-field', '--no-plots',
        str(fits_path),
        '--uniformize', '0',
        '--cpulimit',   '60',
        '--ra',         str(ra_dd),
        '--dec',        str(dec_dd),
        '--radius',     str(args.radius),
        '--new-fits',   str(outfile),
        '--overwrite',
        '--downsample', str(args.downsample),
    ]

    if args.small_field:
        cmd_base += ['-L', '1.5', '-H', '5', '-u', 'amw']

    if args.tweak_order != 0:
        tweak_flags = ['--tweak-order', str(args.tweak_order)]
    else:
        tweak_flags = ['-T']

    _run_astrometry(cmd_base, tweak_flags, fits_path, log_astro)

    log_astro.close()

    if not outfile.exists():
        print(bcolors.FAIL
              + f'ERROR: astrometry.net failed on {fits_path}'
              + bcolors.ENDC)
        log_astro.name  # already closed, just for reference
        return

    print(bcolors.OKGREEN
          + f'WCS calibration successful: {outfile}'
          + bcolors.ENDC)

    # Clean up astrometry.net temporary files
    cwd = Path.cwd()
    workdir = fits_path.parent if fits_path.parent != Path('.') else cwd

    for pattern in ('*.axy', '*indx.xyls', '*.corr', '*.match',
                    '*.rdls', '*.solved', '*.wcs'):
        for f in workdir.glob(pattern):
            try:
                f.unlink()
            except OSError:
                pass

    # Convert SIP distortion keywords to PV format
    ok = sip_to_pv.sip_to_pv(
        infile     = str(outfile),
        outfile    = str(outfile),
        tpv_format = True)

    if not ok:
        print(bcolors.WARNING
              + f'sip_to_pv failed. {outfile} retains SIP keywords only.'
              + bcolors.ENDC)

    log_sip2pv.close()


def _run_astrometry(
    cmd_base: list[str],
    tweak_flags: list[str],
    fits_path: Path,
    log_file,
) -> int:
    """Run solve-field using astrometry.net's built-in image2xy source detector.

    Parameters
    ----------
    cmd_base : list[str]
        Base solve-field command (without tweak flags).
    tweak_flags : list[str]
        ``['--tweak-order', 'N']`` or ``['-T']``.
    fits_path : Path
        Input FITS path — used only for the success message.
    log_file : file
        Open log file; stderr from solve-field is appended on failure.

    Returns
    -------
    int
        Process return code (0 = success).
    """
    cmd = cmd_base + tweak_flags
    print(' '.join(cmd))

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        print(bcolors.OKGREEN
              + f'WCS calibration successful: {fits_path}'
              + bcolors.ENDC)
    else:
        log_file.write(proc.stderr)

    return proc.returncode


if __name__ == '__main__':
    main(sys.argv[1:])
