"""
fits_tools.py — FITS file utilities and WCS coordinate transformations.

Provides coordinate conversion (HMS ↔ decimal degrees), batch sky-to-pixel
transformations, and pixel-scale extraction from WCS headers.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"


import sys
from pathlib import Path

import numpy as np
from astropy import coordinates as coord, units as u, wcs
from astropy.io import ascii, fits
from astropy.wcs.utils import proj_plane_pixel_scales

from utils import bcolors


def convert_hms_dd(
    ra: str | float,
    dec: str | float,
) -> tuple[float, float]:
    """Convert RA/Dec from sexagesimal strings or decimal floats to decimal degrees.

    Parameters
    ----------
    ra : str | float
        Right ascension.  Sexagesimal if it contains ``':'``.
    dec : str | float
        Declination.  Sexagesimal if it contains ``':'``.

    Returns
    -------
    tuple[float, float]
        ``(ra_deg, dec_deg)`` in decimal degrees.
    """
    ra_str  = str(ra).strip().strip('"')
    dec_str = str(dec).strip().strip('"')

    if ':' in ra_str or ':' in dec_str:
        try:
            c = coord.SkyCoord(ra_str, dec_str,
                               unit=(u.hour, u.degree), frame='icrs')
            return float(c.ra.deg), float(c.dec.deg)
        except Exception as exc:
            print(bcolors.FAIL
                  + f'Cannot parse coordinates {ra_str!r} {dec_str!r}: {exc}'
                  + bcolors.ENDC)
            sys.exit(1)

    try:
        return float(ra_str), float(dec_str)
    except ValueError as exc:
        print(bcolors.FAIL
              + f'Cannot parse coordinates {ra_str!r} {dec_str!r}: {exc}'
              + bcolors.ENDC)
        sys.exit(1)


def get_header(
    filepath: str | Path,
    keyword: str,
) -> object:
    """Read a single keyword value from the primary FITS header.

    Parameters
    ----------
    filepath : str | Path
        Path to the FITS file.
    keyword : str
        Header keyword.

    Returns
    -------
    object
        Header value (type depends on the keyword).
    """
    return fits.getheader(str(filepath))[keyword]


def pix2arcsec(fits_path: str | Path) -> float:
    """Return the median pixel scale of a FITS image in arcsec/pixel.

    Parameters
    ----------
    fits_path : str | Path
        Path to the FITS file.

    Returns
    -------
    float
        Median pixel scale in arcsec/pixel.
    """
    with fits.open(str(fits_path)) as hdu:
        header = hdu[0].header.copy()
        if len(hdu) > 1:
            header += hdu[1].header

    return float(np.median(proj_plane_pixel_scales(wcs.WCS(header))) * 3600.0)


def _wcs_from_fits(
    fits_path: str | Path,
) -> tuple[wcs.WCS, int, int]:
    """Load a 2-D WCS and image dimensions from a FITS file.

    Parameters
    ----------
    fits_path : str | Path
        Path to the FITS file.

    Returns
    -------
    tuple[wcs.WCS, int, int]
        ``(wcs_obj, naxis1, naxis2)``.
    """
    with fits.open(str(fits_path)) as hdulist:
        for hdu in hdulist:
            try:
                w = wcs.WCS(hdu.header)
                if w.array_shape is not None and len(w.array_shape) == 2:
                    naxis2, naxis1 = w.array_shape
                    return w, naxis1, naxis2
            except Exception:
                continue

    raise ValueError(f'No valid 2-D WCS found in {fits_path}')


def sky2xy(
    fits_path: str | Path,
    ra: float | None = None,
    dec: float | None = None,
    catalog: str | Path | None = None,
) -> tuple[float, float] | np.ndarray | bool:
    """Convert sky coordinates to 1-indexed pixel positions.

    Single-object mode (pass *ra* and *dec*):
        Returns ``(x, y)`` or ``False`` if the position is outside the image.

    Multi-object mode (pass *catalog*):
        Reads all positions from a two-column ASCII file and applies a
        vectorised batch WCS transformation.  Positions outside the image
        footprint are silently excluded.

    All returned pixel positions follow the 1-indexed FITS convention.

    Parameters
    ----------
    fits_path : str | Path
        Path to the FITS file.
    ra : float | None
        Right ascension in decimal degrees (single-object mode).
    dec : float | None
        Declination in decimal degrees (single-object mode).
    catalog : str | Path | None
        Path to a two-column (RA, Dec) ASCII file (multi-object mode).

    Returns
    -------
    tuple[float, float] | np.ndarray | bool
        ``(x, y)`` in single-object mode, ``np.ndarray`` of shape ``(N, 2)``
        in multi-object mode, or ``False`` if outside the footprint (single mode).
    """
    w, naxis1, naxis2 = _wcs_from_fits(fits_path)

    if catalog is None:
        x, y = w.wcs_world2pix(float(ra), float(dec), 1)
        if np.isnan(x) or np.isnan(y):
            x, y = w.wcs_world2pix(float(dec), float(ra), 1)
        if not (1 <= x <= naxis1 and 1 <= y <= naxis2):
            return False
        return float(x), float(y)

    # Vectorised batch transformation
    cat     = ascii.read(str(catalog))
    keys    = cat.colnames
    ra_arr  = np.asarray(cat[keys[0]], dtype=float)
    dec_arr = np.asarray(cat[keys[1]], dtype=float)

    x_arr, y_arr = w.wcs_world2pix(ra_arr, dec_arr, 1)
    in_frame     = ((x_arr >= 1) & (x_arr <= naxis1)
                    & (y_arr >= 1) & (y_arr <= naxis2)
                    & np.isfinite(x_arr) & np.isfinite(y_arr))

    return np.column_stack([x_arr[in_frame], y_arr[in_frame]])


def xy2sky(
    fits_path: str | Path,
    x: float,
    y: float,
    sep: str = ':',
) -> tuple[str, str]:
    """Convert 1-indexed pixel coordinates to sexagesimal sky coordinates.

    Parameters
    ----------
    fits_path : str | Path
        Path to the FITS file.
    x : float
        x pixel position (1-indexed).
    y : float
        y pixel position (1-indexed).
    sep : str
        Delimiter for HMS/DMS output (default ``':'``).

    Returns
    -------
    tuple[str, str]
        ``(ra_string, dec_string)`` in HMS and DMS format.
    """
    w       = wcs.WCS(fits.getheader(str(fits_path)))
    sky_deg = w.wcs_pix2world([[float(x), float(y)]], 1)[0] * u.deg
    c       = coord.SkyCoord(sky_deg[0], sky_deg[1], frame='fk5')
    parts   = c.to_string(style='hmsdms', sep=sep, precision=3).split()
    return parts[0], parts[1]
