#!/usr/bin/env python
"""
improve_astrometry.py — Improve the WCS calibration of a FITS image.

Uses astrometry.net's ``solve-field`` for blind astrometric calibration.
Source positions are detected with ``sep`` (Python/C, no external binary)
and written to a FITS BINTABLE (xylist).  This is passed to ``solve-field``
via ``--xylist``, which bypasses all Python helper scripts inside the
astrometry.net installation (``image2pnm``, ``removelines``) that require
a Python 3.8 build of the ``astrometry`` package and cannot be used when
the pipeline runs on Python 3.11.

After a successful solve, the WCS keywords from the ``.wcs`` output file
are merged into the original FITS image and written as ``<stem>_wcs.fits``.
The SIP distortion terms are then converted to PV format with sip_to_pv.

Usage
-----
    python improve_astrometry.py --ra 11:33:41.550 --dec 00:43:33.50 \\
        --fits SN2015bn_SDSS_r_wcs.fits
"""

from __future__ import annotations

__version__ = "2026-05-30"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import   argparse
import   subprocess
import   sys
from     pathlib import Path

import   numpy as np
import   sep
from     astropy.io import fits
from     astropy import wcs as astropy_wcs

import   fits_tools
import   sip_to_pv
from     utils import bcolors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_parser() -> argparse.ArgumentParser:
    """Build the argument parser for improve_astrometry."""
    p = argparse.ArgumentParser(
        description=(
            'Improve WCS calibration with astrometry.net '
            'and convert SIP distortion terms to PV format.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--fits',        type=str, required=True,
                   help='Input FITS filename.')
    p.add_argument('--ra',          type=str, required=True,
                   help='RA(J2000) in HMS or decimal degrees.')
    p.add_argument('--dec',         type=str, required=True,
                   help='Dec(J2000) in HMS or decimal degrees.')
    p.add_argument('--radius',      type=float, default=0.125,
                   help='Search radius around field centre (deg).')
    p.add_argument('--det-thresh',  type=float, default=5.0,
                   help='Source detection threshold (sigma) for xylist generation.')
    p.add_argument('--max-sources', type=int,   default=500,
                   help='Maximum number of sources passed to solve-field.')
    p.add_argument('--small-field', action='store_true', default=False,
                   help='Set for fields smaller than 5 arcmin.')
    p.add_argument('--tweak-order', type=int,   default=2,
                   help='Polynomial order for SIP WCS corrections (0 = no SIP).')
    return p


# ---------------------------------------------------------------------------
# Source detection → xylist
# ---------------------------------------------------------------------------

def _make_xylist(
    fits_path: Path,
    detect_thresh: float = 5.0,
    max_sources: int = 500,
) -> tuple[Path, int, int]:
    """Detect sources with sep and write their positions to a FITS xylist.

    astrometry.net's ``solve-field`` accepts a pre-computed source list via
    ``--xylist``, which avoids running ``image2pnm`` / ``removelines`` — the
    helper scripts that require a Python 3.8 ``astrometry`` package and crash
    on Python 3.11.

    Parameters
    ----------
    fits_path : Path
        Input science FITS image.
    detect_thresh : float
        Detection threshold in sigma above background (default: 5.0).
    max_sources : int
        Maximum number of sources to include (brightest-first; default: 500).

    Returns
    -------
    tuple[Path, int, int]
        ``(xylist_path, image_width, image_height)`` — path to the xylist
        FITS file and image dimensions needed for ``--width``/``--height``.
    """
    from astropy.table import Table

    with fits.open(str(fits_path)) as hdulist:
        data = hdulist[0].data
        if data is None and len(hdulist) > 1:
            data = hdulist[1].data

    data = np.ascontiguousarray(data.astype(np.float64))
    naxis2, naxis1 = data.shape

    # Background estimation
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    data_sub = np.ascontiguousarray(data - bkg)

    # Source extraction
    thresh = detect_thresh * bkg.globalrms
    sep.set_extract_pixstack(min(2_000_000, max(500_000, data_sub.size // 3)))
    sep.set_sub_object_limit(262144)

    objects = sep.extract(data_sub, thresh, minarea=5, filter_type='matched')

    if len(objects) == 0:
        raise RuntimeError(
            f'No sources detected at {detect_thresh}σ in {fits_path}. '
            'Lower --det-thresh or check the image.')

    # Sort brightest-first, limit to max_sources
    idx = np.argsort(objects['flux'])[::-1][:max_sources]

    # FITS BINTABLE with 1-indexed X, Y positions
    xylist_path = fits_path.with_suffix('.xyls')
    Table({
        'X':    objects['x'][idx] + 1.0,   # sep 0-indexed → FITS 1-indexed
        'Y':    objects['y'][idx] + 1.0,
        'FLUX': objects['flux'][idx],
    }).write(str(xylist_path), format='fits', overwrite=True)

    print(bcolors.OKGREEN
          + f'{len(idx)} sources written to {xylist_path.name}'
          + bcolors.ENDC)

    return xylist_path, naxis1, naxis2


# ---------------------------------------------------------------------------
# Apply WCS solution to original image
# ---------------------------------------------------------------------------

def _apply_wcs(original_fits: Path, wcs_file: Path, output_fits: Path) -> None:
    """Merge the solved WCS from wcs_file into original_fits and save.

    Strips all existing WCS keywords from the original header, then copies
    the full WCS solution (including SIP distortion terms) from the
    astrometry.net output ``.wcs`` file.

    Parameters
    ----------
    original_fits : Path
        Original science FITS image (source of pixel data and non-WCS header).
    wcs_file : Path
        WCS FITS file produced by solve-field (contains only WCS keywords).
    output_fits : Path
        Output path for the WCS-corrected FITS image.
    """
    # WCS keyword prefixes to strip from the original header
    _WCS_PREFIXES = (
        'WCSAXES', 'CTYPE', 'CRVAL', 'CRPIX', 'CD1_', 'CD2_',
        'CDELT', 'CROTA', 'LONPOLE', 'LATPOLE', 'EQUINOX', 'RADESYS',
        'A_', 'B_', 'AP_', 'BP_',   # SIP coefficients
        'PV',                         # PV distortion
        'PC',                         # PC matrix
    )

    with fits.open(str(original_fits)) as orig:
        img_data   = orig[0].data.copy()
        img_header = orig[0].header.copy()

    # Remove stale WCS keywords
    for key in list(img_header.keys()):
        if any(key.startswith(pfx) for pfx in _WCS_PREFIXES):
            del img_header[key]

    # Merge solved WCS
    wcs_header = fits.getheader(str(wcs_file))
    for key, val, *comment in wcs_header.cards:
        if key in ('', 'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'EXTEND', 'END'):
            continue
        img_header[key] = (val, comment[0] if comment else '')

    fits.writeto(str(output_fits), img_data, img_header,
                 overwrite=True, output_verify='silentfix')


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(args: list[str] | None = None) -> None:
    """Run the astrometry improvement pipeline.

    Parameters
    ----------
    args : list[str] | None
        Command-line arguments. ``sys.argv[1:]`` when None.
    """
    parser = get_parser()
    args   = parser.parse_args(args)

    fits_path = Path(args.fits)
    outfile   = fits_path.with_name(fits_path.stem + '_wcs' + fits_path.suffix)

    log_astro  = open('astro.log',       'w')
    log_sip2pv = open('astro_sip2pv.log', 'w')

    ra_dd, dec_dd = fits_tools.convert_hms_dd(args.ra, args.dec)

    print(bcolors.HEADER + 'Process image: ' + str(fits_path) + bcolors.ENDC)

    # Sanitise FITS headers before processing
    with fits.open(str(fits_path)) as hdulist:
        hdu_data   = hdulist[0].data
        hdu_header = hdulist[0].header
    fits.writeto(str(fits_path), hdu_data, hdu_header,
                 overwrite=True, output_verify='silentfix')

    # 1. Detect sources with sep and write xylist
    try:
        xylist, img_w, img_h = _make_xylist(
            fits_path,
            detect_thresh = args.det_thresh,
            max_sources   = args.max_sources)
    except RuntimeError as exc:
        print(bcolors.FAIL + str(exc) + bcolors.ENDC)
        return

    # 2. Build solve-field command using the xylist
    #    --xylist bypasses image2pnm / removelines (Python 3.8 helpers that
    #    crash on Python 3.11) so solve-field runs as a pure C binary.
    # Pass the xylist as a positional argument (solve-field auto-detects FITS
    # BINTABLE vs image from the file content).
    # -w / --width and -e / --height tell solve-field the image dimensions.
    cmd = [
        'solve-field', '--no-plots',
        '--no-remove-lines',    # skip removelines (Python 3.8 helper, breaks on 3.11)
        str(xylist),
        '-w', str(img_w),
        '-e', str(img_h),
        '--uniformize', '0',
        '--cpulimit',   '60',
        '--ra',         str(ra_dd),
        '--dec',        str(dec_dd),
        '--radius',     str(args.radius),
        '--overwrite',
    ]

    if args.small_field:
        cmd += ['-L', '1.5', '-H', '5', '-u', 'amw']

    if args.tweak_order != 0:
        cmd += ['--tweak-order', str(args.tweak_order)]
    else:
        cmd += ['-T']

    print(' '.join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        log_astro.write(proc.stderr)
        log_astro.close()
        print(bcolors.FAIL
              + f'ERROR: astrometry.net failed on {fits_path}'
              + bcolors.ENDC)
        return

    print(bcolors.OKGREEN + 'Astrometry solved.' + bcolors.ENDC)
    log_astro.close()

    # 3. Apply WCS solution to original image
    #    solve-field writes the WCS to <xylist_stem>.wcs
    wcs_file = xylist.with_suffix('.wcs')
    if not wcs_file.exists():
        print(bcolors.FAIL
              + f'ERROR: WCS output file not found: {wcs_file}'
              + bcolors.ENDC)
        return

    _apply_wcs(fits_path, wcs_file, outfile)
    print(bcolors.OKGREEN + f'WCS applied → {outfile}' + bcolors.ENDC)

    # 4. Clean up astrometry.net temporary files
    workdir = fits_path.parent if str(fits_path.parent) != '.' else Path.cwd()
    for pattern in ('*.axy', '*indx.xyls', '*.corr', '*.match',
                    '*.rdls', '*.solved', '*.wcs', '*.xyls'):
        for f in workdir.glob(pattern):
            try:
                f.unlink()
            except OSError:
                pass

    # 5. Convert SIP distortion keywords to PV format
    ok = sip_to_pv.sip_to_pv(
        infile     = str(outfile),
        outfile    = str(outfile),
        tpv_format = True)

    if not ok:
        print(bcolors.WARNING
              + f'sip_to_pv failed. {outfile} retains SIP keywords only.'
              + bcolors.ENDC)

    log_sip2pv.close()


if __name__ == '__main__':
    main(sys.argv[1:])
