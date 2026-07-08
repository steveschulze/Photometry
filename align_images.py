#!/usr/bin/env python
"""
align_images.py — Align the WCS astrometry of a science image to a reference.

Uses ``sep`` for source detection and ``SCAMP`` for astrometric alignment.
SCAMP matches sources detected in the new image against a reference catalog
built from the reference image and updates the WCS keywords of the new image.

Source positions use windowed centroids (equivalent to SExtractor's
``XWIN_IMAGE`` / ``YWIN_IMAGE``) for sub-pixel accuracy.

The aligned image is saved as ``<new-image-stem>_astro.fits``.
Diagnostic quality metrics from SCAMP are written to
``<new-image-stem>_alignment.log``.

Usage
-----
    python align_images.py \\
        --ref-image reference.fits --new-image science.fits
"""

from __future__ import annotations

__version__ = "2026-06-01"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import   argparse
import   logging
import   subprocess
import   sys
from     pathlib import Path

import   numpy as np
import   sep
from     astropy.io import fits, votable
from     astropy import wcs as astropy_wcs

from     utils import bcolors, setup_logging


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_parser() -> argparse.ArgumentParser:
    """Build the argument parser for align_images."""
    p = argparse.ArgumentParser(
        description='Align WCS astrometry of a science image to a reference.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--ref-image',
                   type=str, required=True,
                   help='Reference image (must have a valid astrometric WCS).')
    p.add_argument('--new-image',
                   type=str, required=True,
                   help='Science image to be aligned.')
    p.add_argument('--det-thresh',
                   type=float, default=5.0,
                   help='Source detection threshold in sigma above background.')
    p.add_argument('--distort-degrees',
                   type=int, default=1,
                   help='SCAMP polynomial degree for the astrometric distortion model.')
    p.add_argument('--position-maxerr',
                   type=float, default=0.15,
                   help='Maximum positional offset (deg) for SCAMP source matching.')
    p.add_argument('--sn-thresholds',
                   type=str, default='10.0,40.0',
                   help='SCAMP S/N thresholds (low,high) for two-pass matching.')
    p.add_argument('--keep-temp',
                   action='store_true', default=False,
                   help='Keep intermediate LDAC catalogs and SCAMP output files.')
    p.add_argument('--loglevel',
                   type=str, default='INFO',
                   help='Logging level (DEBUG/INFO/WARNING/ERROR/CRITICAL).')
    return p


# ---------------------------------------------------------------------------
# LDAC catalog creation
# ---------------------------------------------------------------------------

def _make_ldac(
    fits_path: Path,
    ldac_path: Path,
    detect_thresh: float = 5.0,
    world_coords: bool = False,
    logger: logging.Logger | None = None,
) -> int:
    """Detect sources with sep and write a FITS LDAC catalog for SCAMP.

    Creates a standard FITS LDAC file with two binary-table extensions:

    * ``LDAC_IMHEAD`` — the original image header serialised as a single
      character column; SCAMP reads this to recover pixel scale and image
      geometry.
    * ``LDAC_OBJECTS`` — source table with positional, morphological, and
      photometric columns.

    Two catalog modes are supported:

    * ``world_coords=True`` — reference catalog: positions as
      ``ALPHA_J2000`` / ``DELTA_J2000`` (degrees).  Requires a valid WCS.
    * ``world_coords=False`` — input catalog: positions as
      ``XWIN_IMAGE`` / ``YWIN_IMAGE`` (1-indexed pixels), matching the
      default ``CENTROID_KEYS`` expected by SCAMP.

    Windowed centroids (``sep.winpos``) are used throughout to improve
    sub-pixel position accuracy, equivalent to SExtractor's ``XWIN_IMAGE`` /
    ``YWIN_IMAGE``.

    Parameters
    ----------
    fits_path : Path
        Input FITS image.
    ldac_path : Path
        Output LDAC FITS file.
    detect_thresh : float
        Detection threshold in sigma above background (default: 5.0).
    world_coords : bool
        When True, write RA/Dec columns for a SCAMP reference catalog.
    logger : logging.Logger | None
        Logger instance for progress and diagnostics.

    Returns
    -------
    int
        Number of sources written to the catalog.

    Raises
    ------
    RuntimeError
        If no sources are detected or (for world_coords=True) the image
        has no valid WCS.
    """
    with fits.open(str(fits_path)) as hdulist:
        img_header = hdulist[0].header.copy()
        data       = hdulist[0].data
        if data is None and len(hdulist) > 1:
            data = hdulist[1].data

    data = np.ascontiguousarray(data.astype(np.float64))

    # Background estimation
    bkg     = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    data_sub = np.ascontiguousarray(data - bkg)

    # Source detection
    thresh = detect_thresh * bkg.globalrms
    sep.set_extract_pixstack(min(2_000_000, max(500_000, data_sub.size // 3)))
    sep.set_sub_object_limit(262144)
    objects = sep.extract(data_sub, thresh, minarea=5, filter_type='matched')

    n = len(objects)
    if n == 0:
        msg = f'No sources detected at {detect_thresh}σ in {fits_path.name}'
        if logger:
            logger.error(msg)
        raise RuntimeError(msg)

    # Clamp theta to [-π/2, π/2].  sep.extract can return values fractionally
    # outside this range due to floating-point drift in second-moment
    # computation; sep.sum_ellipse raises "invalid aperture parameters" for
    # even one ULP beyond the boundary.
    objects['theta'] = np.clip(objects['theta'], -np.pi / 2.0, np.pi / 2.0)

    if logger:
        logger.info('  %s: %d sources at %.1fσ (threshold %.4f ADU)',
                    fits_path.name, n, detect_thresh, thresh)

    # Windowed centroids (XWIN / YWIN equivalent)
    sig = np.maximum(objects['a'] / np.sqrt(2.0), 0.5)
    try:
        xwin, ywin, flag_win = sep.winpos(data_sub, objects['x'], objects['y'], sig)
        good_win = flag_win == 0
        x_pos    = np.where(good_win, xwin, objects['x'])
        y_pos    = np.where(good_win, ywin, objects['y'])
        if logger:
            logger.info('  Windowed centroids: %d/%d succeeded', int(good_win.sum()), n)
    except Exception as exc:
        if logger:
            logger.warning('sep.winpos failed (%s); using isophotal centroids', exc)
        x_pos = objects['x']
        y_pos = objects['y']

    # SNR estimate: flux / (background RMS × sqrt(source area in pixels)).
    # sep.extract() only populates 'fluxerr' when an error array is passed;
    # without one (absolute-threshold mode) we derive the noise from the
    # background and the segmentation footprint size.
    flux_err   = bkg.globalrms * np.sqrt(np.maximum(objects['npix'].astype(float), 1.0))
    snr        = np.maximum(objects['flux'] / np.maximum(flux_err, 1e-10), 1.0)
    pos_err_px = np.clip(objects['a'] / snr, 0.05, 10.0).astype(np.float32)

    # ------------------------------------------------------------------
    # LDAC_IMHEAD — image header as a single long character string
    # ------------------------------------------------------------------
    header_str = img_header.tostring()
    # Pad to a multiple of 2880 bytes (FITS block size) as SCAMP expects
    block = 2880
    pad   = (block - len(header_str) % block) % block
    header_str += ' ' * pad

    hdr_col = fits.Column(
        name='Field Header Card',
        format=f'{len(header_str)}A',
        array=[header_str])
    ldac_head_hdu      = fits.BinTableHDU.from_columns([hdr_col])
    ldac_head_hdu.name = 'LDAC_IMHEAD'

    # ------------------------------------------------------------------
    # LDAC_OBJECTS — source table
    # ------------------------------------------------------------------
    if world_coords:
        # Reference catalog: RA/Dec in degrees
        try:
            w   = astropy_wcs.WCS(img_header)
            sky = w.pixel_to_world(x_pos, y_pos)
            ra  = np.asarray(sky.ra.deg,  dtype=np.float64)
            dec = np.asarray(sky.dec.deg, dtype=np.float64)
        except Exception as exc:
            raise RuntimeError(
                f'WCS pixel→world conversion failed for {fits_path.name}: {exc}'
            ) from exc

        # Positional errors in degrees (pixel error × pixel scale)
        try:
            from astropy.wcs.utils import proj_plane_pixel_scales
            pix_scale = float(np.median(proj_plane_pixel_scales(w)) * 3600.0)
        except Exception:
            pix_scale = 1.0

        pos_err_deg = (pos_err_px * pix_scale / 3600.0).astype(np.float32)

        # Instrumental magnitudes for SCAMP S/N selection
        flux_aper, _, _ = sep.sum_circle(data_sub, x_pos, y_pos, 5.0)
        good_flux = flux_aper > 0
        mag_auto  = np.full(n, 99.0, dtype=np.float32)
        mag_err   = np.full(n, 99.0, dtype=np.float32)
        mag_auto[good_flux] = (-2.5 * np.log10(flux_aper[good_flux])).astype(np.float32)
        mag_err[good_flux]  = (2.5 / np.log(10) / snr[good_flux]).astype(np.float32)

        cols = [
            fits.Column(name='X_WORLD',        format='1D', array=ra),
            fits.Column(name='Y_WORLD',        format='1D', array=dec),
            fits.Column(name='ERRA_WORLD',     format='1E', array=pos_err_deg),
            fits.Column(name='ERRB_WORLD',     format='1E', array=pos_err_deg),
            fits.Column(name='ERRTHETA_WORLD', format='1E', array=np.zeros(n, np.float32)),
            fits.Column(name='MAG_AUTO',       format='1E', array=mag_auto),
            fits.Column(name='MAGERR_AUTO',    format='1E', array=mag_err),
        ]
    else:
        # Input catalog: 1-indexed pixel positions (SCAMP convention)
        cols = [
            fits.Column(name='XWIN_IMAGE',        format='1E',
                        array=(x_pos + 1.0).astype(np.float32)),
            fits.Column(name='YWIN_IMAGE',        format='1E',
                        array=(y_pos + 1.0).astype(np.float32)),
            fits.Column(name='ERRAWIN_IMAGE',     format='1E', array=pos_err_px),
            fits.Column(name='ERRBWIN_IMAGE',     format='1E', array=pos_err_px),
            fits.Column(name='ERRTHETAWIN_IMAGE', format='1E',
                        array=np.zeros(n, np.float32)),
            fits.Column(name='FLUX_AUTO',         format='1E',
                        array=objects['flux'].astype(np.float32)),
            fits.Column(name='FLUXERR_AUTO',      format='1E',
                        array=flux_err.astype(np.float32)),
            fits.Column(name='FLUX_RADIUS',       format='1E',
                        array=objects['a'].astype(np.float32)),
            fits.Column(name='FLAGS',             format='1J',
                        array=objects['flag'].astype(np.int32)),
        ]

    ldac_obj_hdu      = fits.BinTableHDU.from_columns(cols)
    ldac_obj_hdu.name = 'LDAC_OBJECTS'

    fits.HDUList([fits.PrimaryHDU(), ldac_head_hdu, ldac_obj_hdu]).writeto(
        str(ldac_path), overwrite=True)

    if logger:
        logger.info('  LDAC written: %s (%d sources)', ldac_path.name, n)

    return n


# ---------------------------------------------------------------------------
# SCAMP quality assessment
# ---------------------------------------------------------------------------

def _parse_scamp_xml(
    xml_path: Path,
    logger: logging.Logger | None = None,
) -> dict:
    """Parse SCAMP's VOTable XML output for astrometric quality metrics.

    Quality indicators extracted from the astrometry statistics table:

    * ``XY_Contrast`` — pattern-matching contrast; values > 3 indicate a
      reliable, unambiguous solution.
    * ``DX``, ``DY`` — mean residual offsets in RA*cos(Dec) and Dec (arcsec).
    * ``Sigma_Alpha_Int``, ``Sigma_Delta_Int`` — internal 1-sigma scatter of
      matched pairs in RA and Dec (arcsec); the primary figure of merit for
      the internal astrometric precision.
    * ``Chi2_Int`` — reduced chi-squared of the internal astrometric fit;
      values close to 1 indicate a well-conditioned solution.
    * ``Nb_Sources``, ``Nb_Matched`` — total and matched source counts.

    A known bug in some SCAMP versions writes ``datatype="*"`` which is
    invalid VOTable syntax.  This is patched transparently before parsing.

    Parameters
    ----------
    xml_path : Path
        SCAMP VOTable XML output file.
    logger : logging.Logger | None
        Logger instance.

    Returns
    -------
    dict
        Parsed quality metrics.  Empty dict if parsing fails or the file
        does not exist.
    """
    if not xml_path.exists():
        if logger:
            logger.warning('SCAMP XML not found: %s', xml_path)
        return {}

    # Patch invalid datatype="*" produced by some SCAMP versions
    try:
        xml_text = xml_path.read_text(encoding='utf-8', errors='replace')
        if 'datatype="*"' in xml_text:
            xml_path.write_text(
                xml_text.replace('datatype="*"', 'datatype="char"'),
                encoding='utf-8')
            if logger:
                logger.debug('Patched invalid datatype="*" in SCAMP XML')
    except OSError:
        pass

    try:
        vot = votable.parse(str(xml_path), invalid='mask')
        tbl = vot.get_first_table().array
    except Exception as exc:
        if logger:
            logger.warning('Failed to parse SCAMP XML: %s', exc)
        return {}

    def _scalar(field: str) -> float | None:
        try:
            v = tbl[field]
            return float(v[0]) if len(v) else None
        except (KeyError, ValueError, IndexError):
            return None

    def _string(field: str) -> str | None:
        try:
            v = tbl[field]
            return str(v[0]).strip() if len(v) else None
        except (KeyError, ValueError, IndexError):
            return None

    # DX / DY are offsets in degrees; convert to arcsec
    dx_deg = _scalar('DX') or 0.0
    dy_deg = _scalar('DY') or 0.0

    metrics: dict = {
        'catalog_name':    _string('Catalog_Name'),
        'nb_sources':      _scalar('Nb_Sources'),
        'nb_matched':      _scalar('Nb_Matched'),
        'xy_contrast':     _scalar('XY_Contrast'),
        'dx_arcsec':       dx_deg * 3600.0,
        'dy_arcsec':       dy_deg * 3600.0,
        'sigma_alpha_int': _scalar('Sigma_Alpha_Int'),
        'sigma_delta_int': _scalar('Sigma_Delta_Int'),
        'chi2_int':        _scalar('Chi2_Int'),
    }

    return {k: v for k, v in metrics.items() if v is not None}


def _log_metrics(metrics: dict, logger: logging.Logger | None = None) -> None:
    """Print and log SCAMP quality metrics in a formatted table.

    Parameters
    ----------
    metrics : dict
        Quality metrics dict as returned by :func:`_parse_scamp_xml`.
    logger : logging.Logger | None
        Logger instance.
    """
    rows = [
        ('nb_sources',      'Sources detected',         '',       '.0f'),
        ('nb_matched',      'Sources matched',          '',       '.0f'),
        ('xy_contrast',     'XY contrast',              '',       '.2f'),
        ('dx_arcsec',       'Mean dRA',                 'arcsec', '.3f'),
        ('dy_arcsec',       'Mean dDec',                'arcsec', '.3f'),
        ('sigma_alpha_int', 'Sigma_alpha (internal)',   'arcsec', '.3f'),
        ('sigma_delta_int', 'Sigma_delta (internal)',   'arcsec', '.3f'),
        ('chi2_int',        'Chi2 (internal)',          '',       '.3f'),
    ]

    for key, label, unit, fmt in rows:
        val = metrics.get(key)
        if val is None:
            continue
        val_str = format(val, fmt)
        line    = f'  {label:<28s} {val_str:>10s}  {unit}'.rstrip()
        print(line)
        if logger:
            logger.info(line.strip())


# ---------------------------------------------------------------------------
# SCAMP configuration file
# ---------------------------------------------------------------------------

def _make_scamp_config(config_path: Path, logger: logging.Logger | None = None) -> bool:
    """Generate a default SCAMP configuration file via ``scamp -dd``.

    SCAMP 2.x looks for a ``scamp.conf`` in the current directory and emits
    a warning if none is found, then falls back to internal defaults that may
    differ from the expected behaviour.  Generating an explicit config file
    ensures reproducible settings and silences the warning.

    Parameters
    ----------
    config_path : Path
        Destination path for the generated configuration file.
    logger : logging.Logger | None
        Logger instance.

    Returns
    -------
    bool
        True if the config file was written, False if ``scamp -dd`` failed
        (in which case SCAMP will use its internal defaults).
    """
    proc = subprocess.run(['scamp', '-dd'], capture_output=True, text=True)
    if proc.returncode != 0 or not proc.stdout.strip():
        if logger:
            logger.warning('scamp -dd failed; SCAMP will use internal defaults')
        return False
    config_path.write_text(proc.stdout, encoding='utf-8')
    if logger:
        logger.debug('SCAMP config written to %s', config_path.name)
    return True


# ---------------------------------------------------------------------------
# Apply SCAMP solution
# ---------------------------------------------------------------------------

def _apply_scamp_solution(
    original_fits: Path,
    head_file: Path,
    output_fits: Path,
    logger: logging.Logger | None = None,
) -> None:
    """Replace the WCS in original_fits with SCAMP's updated solution.

    Strips all existing WCS keywords from the original header, then applies
    the updated WCS written by SCAMP to the ``.head`` output file.

    SCAMP writes ``.head`` files in FITS-card format: 80-character records,
    potentially separated by newlines depending on the SCAMP version.  The
    parser normalises both styles.

    Parameters
    ----------
    original_fits : Path
        Original science FITS image.
    head_file : Path
        SCAMP ``.head`` file containing the updated WCS keywords.
    output_fits : Path
        Destination path for the aligned FITS image.
    logger : logging.Logger | None
        Logger instance.
    """
    # Keywords to strip from the original header before applying the solution
    _WCS_PREFIXES = (
        'WCSAXES', 'CTYPE', 'CRVAL', 'CRPIX',
        'CD1_', 'CD2_', 'CDELT', 'CROTA',
        'LONPOLE', 'LATPOLE', 'EQUINOX', 'RADESYS',
        'A_', 'B_', 'AP_', 'BP_',
        'PV', 'PC',
    )
    # Whitelist of prefixes to *apply* from the SCAMP .head file.
    # Using a whitelist (rather than a skip-list) avoids accidentally copying
    # non-WCS metadata that SCAMP writes — e.g. ORIGIN, AUTHOR — which may
    # contain non-ASCII characters and cause astropy to raise ValueError.
    _WCS_APPLY_PREFIXES = (
        'WCSAXES', 'CTYPE', 'CRVAL', 'CRPIX', 'CUNIT',
        'CD1_', 'CD2_', 'CDELT', 'CROTA',
        'LONPOLE', 'LATPOLE', 'EQUINOX', 'RADESYS', 'RADECSYS',
        'PV', 'PC',
        'FLXSCALE',   # photometric scale factor (when SOLVE_PHOTOM Y)
        'ASTIRMS',    # SCAMP internal scatter keywords
        'ASTRRMS',    # SCAMP reference scatter keywords
    )

    with fits.open(str(original_fits)) as orig:
        img_data   = orig[0].data.copy()
        img_header = orig[0].header.copy()

    # Remove stale WCS keywords from original header
    for key in list(img_header.keys()):
        if any(key.startswith(pfx) for pfx in _WCS_PREFIXES):
            del img_header[key]

    # Parse SCAMP .head file: 80-char FITS-card records, strip newlines
    head_content = head_file.read_text(encoding='utf-8', errors='replace')
    flat         = head_content.replace('\n', '')
    n_cards      = len(flat) // 80
    applied      = 0

    for i in range(n_cards):
        card = flat[i * 80:(i + 1) * 80]
        kw   = card[:8].strip()
        # Only apply recognised WCS keywords; skip all other SCAMP metadata
        if not kw or '=' not in card:
            continue
        if not any(kw.startswith(pfx) for pfx in _WCS_APPLY_PREFIXES):
            continue
        val_str = card[9:].split('/')[0].strip().strip("'")
        try:
            val: float | str = float(val_str)
        except ValueError:
            val = val_str
        img_header[kw] = val
        applied += 1

    fits.writeto(str(output_fits), img_data, img_header,
                 overwrite=True, output_verify='silentfix')

    if logger:
        logger.info('SCAMP solution applied (%d keywords) -> %s',
                    applied, output_fits.name)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(args: list[str] | None = None) -> None:
    """Run the image alignment pipeline.

    Parameters
    ----------
    args : list[str] | None
        Command-line arguments; defaults to ``sys.argv[1:]``.
    """
    parser = get_parser()
    args   = parser.parse_args(args)

    ref_path = Path(args.ref_image)
    new_path = Path(args.new_image)
    out_path = new_path.with_name(new_path.stem + '_astro' + new_path.suffix)

    # Intermediate files live alongside the new image
    ref_ldac = new_path.with_name(ref_path.stem + '_ref.ldac')
    new_ldac = new_path.with_name(new_path.stem + '_new.ldac')
    xml_path   = new_path.with_name(new_path.stem + '_scamp.xml')
    scamp_conf = new_path.with_name(new_path.stem + '_scamp.conf')
    head_file  = new_ldac.with_suffix('.head')

    # ── Logging ──────────────────────────────────────────────────────────
    log_path = new_path.with_name(new_path.stem + '_alignment.log')
    logger   = setup_logging('align_images', log_path, level=args.loglevel)
    logger.info('align_images.py')
    logger.info('Command: %s', ' '.join(sys.argv))
    logger.info('Reference image  : %s', ref_path)
    logger.info('New image        : %s', new_path)
    logger.info('Output           : %s', out_path)
    logger.info('det_thresh=%.1f  distort_degrees=%d  position_maxerr=%.3f deg',
                args.det_thresh, args.distort_degrees, args.position_maxerr)

    print(bcolors.HEADER + 'align_images — ' + str(new_path) + bcolors.ENDC)

    # ── Step 1: Build LDAC catalogs ───────────────────────────────────────
    logger.info('Step 1: source detection and LDAC catalog generation')

    try:
        n_ref = _make_ldac(ref_path, ref_ldac,
                           detect_thresh=args.det_thresh,
                           world_coords=True,
                           logger=logger)
        print(bcolors.OKGREEN
              + f'Reference catalog: {n_ref} sources -> {ref_ldac.name}'
              + bcolors.ENDC)
    except RuntimeError as exc:
        print(bcolors.FAIL + 'ERROR: ' + str(exc) + bcolors.ENDC)
        logger.error(str(exc))
        return

    try:
        n_new = _make_ldac(new_path, new_ldac,
                           detect_thresh=args.det_thresh,
                           world_coords=False,
                           logger=logger)
        print(bcolors.OKGREEN
              + f'New-image catalog: {n_new} sources -> {new_ldac.name}'
              + bcolors.ENDC)
    except RuntimeError as exc:
        print(bcolors.FAIL + 'ERROR: ' + str(exc) + bcolors.ENDC)
        logger.error(str(exc))
        return

    # ── Step 2: Run SCAMP ────────────────────────────────────────────────
    logger.info('Step 2: astrometric alignment (SCAMP)')

    # Generate a default SCAMP config file so SCAMP doesn't warn about a
    # missing scamp.conf and uses well-defined defaults.
    _make_scamp_config(scamp_conf, logger=logger)

    # Notes on omitted parameters:
    #   -ASTREFRA_KEY / -ASTREFDEC_KEY  — removed; SCAMP 2.x does not
    #     recognise them as valid command-line overrides and warns
    #     "keyword unknown".  X_WORLD/Y_WORLD are the built-in defaults.
    #   -CHECKPLOT_TYPE NONE            — disables checkplot file creation,
    #     preventing SCAMP from calling ImageMagick's `convert` (which may
    #     not be installed and causes a SIGSEGV on some systems).
    cmd = [
        'scamp', str(new_ldac),
        '-c',                str(scamp_conf),
        '-ASTREF_CATALOG',   'FILE',
        '-ASTREFCAT_NAME',   str(ref_ldac),
        '-ASTREFMAG_KEY',    'MAG_AUTO',
        '-CENTROID_KEYS',    'XWIN_IMAGE,YWIN_IMAGE',
        '-CENTROIDERR_KEYS', 'ERRAWIN_IMAGE,ERRBWIN_IMAGE,ERRTHETAWIN_IMAGE',
        '-MATCH',            'Y',
        '-WRITE_XML',        'Y',
        '-XML_NAME',         str(xml_path),
        '-CHECKPLOT_TYPE',   'NONE',
        '-CHECKPLOT_DEV',    'NULL',
        '-DISTORT_DEGREES',  str(args.distort_degrees),
        '-SOLVE_ASTROM',     'Y',
        '-SOLVE_PHOTOM',     'N',
        '-POSITION_MAXERR',  str(args.position_maxerr),
        '-SN_THRESHOLDS',    args.sn_thresholds,
    ]

    logger.info('SCAMP command: %s', ' '.join(cmd))
    print(' '.join(cmd))

    # Let SCAMP output flow to the terminal directly so the user can monitor
    # progress and diagnose matching problems.  stderr is captured so failures
    # are written to the log file.
    proc = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        msg = f'SCAMP failed (exit code {proc.returncode})'
        print(bcolors.FAIL + 'ERROR: ' + msg + bcolors.ENDC)
        logger.error(msg)
        if proc.stderr:
            logger.error('SCAMP stderr:\n%s', proc.stderr.strip())
            print(proc.stderr)
        return

    print(bcolors.OKGREEN + 'SCAMP alignment completed.' + bcolors.ENDC)
    logger.info('SCAMP completed successfully')

    # ── Step 3: Quality assessment ────────────────────────────────────────
    logger.info('Step 3: alignment quality assessment')
    metrics = _parse_scamp_xml(xml_path, logger=logger)

    if metrics:
        print(bcolors.OKBLUE + '\nAlignment quality metrics' + bcolors.ENDC)
        logger.info('=== Alignment quality metrics ===')
        _log_metrics(metrics, logger=logger)
    else:
        logger.warning('No quality metrics extracted from SCAMP XML')

    # ── Step 4: Apply SCAMP solution ──────────────────────────────────────
    logger.info('Step 4: applying SCAMP WCS solution')

    if not head_file.exists():
        msg = f'SCAMP .head file not found: {head_file}'
        print(bcolors.FAIL + 'ERROR: ' + msg + bcolors.ENDC)
        logger.error(msg)
        return

    _apply_scamp_solution(new_path, head_file, out_path, logger=logger)
    print(bcolors.OKGREEN + f'Aligned image written -> {out_path}' + bcolors.ENDC)

    # ── Step 5: Clean up temporary files ─────────────────────────────────
    if not args.keep_temp:
        logger.info('Step 5: removing temporary files')
        for tmp in (ref_ldac, new_ldac, head_file, xml_path, scamp_conf):
            try:
                tmp.unlink(missing_ok=True)
            except OSError:
                pass
    else:
        logger.info('Temporary files kept (--keep-temp)')

    logger.info('Done — output: %s', out_path)


if __name__ == '__main__':
    main(sys.argv[1:])
