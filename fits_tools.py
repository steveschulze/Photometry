"""
fits_tools.py — FITS file utilities and WCS coordinate transformations.

Provides coordinate conversion (HMS ↔ decimal degrees), sky-to-pixel
transformations, and pixel-scale extraction.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

from   astropy import coordinates as coord
from   astropy import wcs
from   astropy.io import ascii, fits
from   astropy import units as u
from   astropy.wcs.utils import proj_plane_pixel_scales
from   misc import bcolors
import numpy as np
import sys


def convert_hms_dd(RA, DEC):
    """Convert RA/Dec from HMS/DMS strings or decimal strings to decimal degrees.

    Accepts either colon-separated sexagesimal strings (``HH:MM:SS.ss``,
    ``±DD:MM:SS.s``) or plain decimal degree strings / floats.

    Parameters
    ----------
    RA : str or float
        Right ascension.  Sexagesimal if it contains ``':'``, otherwise
        interpreted as decimal degrees.
    DEC : str or float
        Declination.  Sexagesimal if it contains ``':'``, otherwise
        interpreted as decimal degrees.

    Returns
    -------
    tuple
        (ra_dd, dec_dd) both in decimal degrees (float).

    Raises
    ------
    SystemExit
        If the coordinate format cannot be parsed.
    """
    ra_str  = str(RA).strip().strip('"')
    dec_str = str(DEC).strip().strip('"')

    has_colon = ':' in ra_str or ':' in dec_str

    if has_colon:
        try:
            c = coord.SkyCoord(ra_str, dec_str,
                               unit=(u.hour, u.degree), frame='icrs')
            return float(c.ra.deg), float(c.dec.deg)
        except Exception as exc:
            print(bcolors.FAIL
                  + f'Cannot parse coordinates "{ra_str}" "{dec_str}": {exc}'
                  + bcolors.ENDC)
            sys.exit(1)

    try:
        return float(ra_str), float(dec_str)
    except ValueError as exc:
        print(bcolors.FAIL
              + f'Cannot parse coordinates "{ra_str}" "{dec_str}": {exc}'
              + bcolors.ENDC)
        sys.exit(1)


def get_header(filepath, keyword):
    """Read a single keyword from a FITS header.

    Parameters
    ----------
    filepath : str
        Path to the FITS file.
    keyword : str
        Header keyword to retrieve.

    Returns
    -------
    value
        Header value (type depends on keyword).
    """
    return fits.getheader(filepath)[keyword]


def pix2arcsec(fits_path):
    """Return the median pixel scale of a FITS image in arcsec/pixel.

    Reads WCS information from the primary (and, if present, first extension)
    header and computes the plate scale using
    :func:`astropy.wcs.utils.proj_plane_pixel_scales`.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.

    Returns
    -------
    float
        Median pixel scale in arcsec/pixel.
    """
    with fits.open(fits_path) as hdu:
        header = hdu[0].header.copy()
        if len(hdu) > 1:
            header += hdu[1].header

    hdu_wcs = wcs.WCS(header)
    return float(np.median(proj_plane_pixel_scales(hdu_wcs)) * 3600.)


def _wcs_from_fits(fits_path):
    """Load a WCS object and image dimensions from a FITS file.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.

    Returns
    -------
    tuple
        (wcs_obj, naxis1, naxis2) — WCS object and image dimensions.
    """
    with fits.open(fits_path) as hdulist:
        # Find the extension that has both WCS and image data
        for hdu in hdulist:
            try:
                w = wcs.WCS(hdu.header)
                if w.array_shape is not None and len(w.array_shape) == 2:
                    naxis2, naxis1 = w.array_shape
                    return w, naxis1, naxis2
            except Exception:
                continue

    raise ValueError(f'No valid 2-D WCS found in {fits_path}')


def sky2xy(FITS, RA=None, DEC=None, CAT=None):
    """Convert sky coordinates to pixel positions for an image.

    When *CAT* is given, reads all positions from a two-column ASCII file
    and uses a **vectorised** WCS transformation (one call for all rows).
    Positions outside the image footprint are silently excluded.

    All returned pixel positions are 1-indexed (FITS convention).

    Parameters
    ----------
    FITS : str
        Path to the FITS file.
    RA : float, optional
        Right ascension in decimal degrees (single-object mode).
    DEC : float, optional
        Declination in decimal degrees (single-object mode).
    CAT : str, optional
        Path to a two-column (RA, Dec) ASCII catalog (multi-object mode).

    Returns
    -------
    tuple or np.ndarray
        Single-object: ``(x, y)`` tuple (1-indexed floats).
        Multi-object: ``np.ndarray`` of shape ``(N, 2)`` with x, y columns.

    Raises
    ------
    SystemExit
        Single-object mode: if object is outside the image footprint or
        WCS transformation fails.
    """
    w, naxis1, naxis2 = _wcs_from_fits(FITS)

    if CAT is None:
        # Single-object mode
        ra_val  = float(RA)
        dec_val = float(DEC)
        x, y = w.wcs_world2pix(ra_val, dec_val, 1)

        if np.isnan(x) or np.isnan(y):
            # Try swapping axes (rare header misconfigurations)
            x, y = w.wcs_world2pix(dec_val, ra_val, 1)

        if x < 1 or x > naxis1 or y < 1 or y > naxis2:
            return False

        return float(x), float(y)

    # Multi-object mode — vectorised batch transformation
    cat_data = ascii.read(CAT)
    keys     = cat_data.colnames

    ra_arr  = np.array(cat_data[keys[0]], dtype=float)
    dec_arr = np.array(cat_data[keys[1]], dtype=float)

    # Batch WCS call: origin=1 → 1-indexed output
    x_arr, y_arr = w.wcs_world2pix(ra_arr, dec_arr, 1)

    # Keep only positions inside image footprint
    in_frame = ((x_arr >= 1) & (x_arr <= naxis1)
                & (y_arr >= 1) & (y_arr <= naxis2)
                & np.isfinite(x_arr) & np.isfinite(y_arr))

    return np.column_stack([x_arr[in_frame], y_arr[in_frame]])


def xy2sky(FITS, X, Y, sep=':'):
    """Convert pixel coordinates to sky coordinates.

    Parameters
    ----------
    FITS : str
        Path to the FITS file.
    X : float
        x pixel position (1-indexed).
    Y : float
        y pixel position (1-indexed).
    sep : str
        Separator for HMS/DMS string output (default ``':'``).

    Returns
    -------
    tuple
        (ra_str, dec_str) formatted as HMS and DMS strings.
    """
    header = fits.getheader(FITS)
    w = wcs.WCS(header)
    sky_deg = w.wcs_pix2world([[float(X), float(Y)]], 1)[0] * u.deg
    c = coord.SkyCoord(sky_deg[0], sky_deg[1], frame='fk5')
    alpha = c.to_string(style='hmsdms', sep=sep, precision=3).split()[0]
    delta = c.to_string(style='hmsdms', sep=sep, precision=2).split()[1]
    return alpha, delta
