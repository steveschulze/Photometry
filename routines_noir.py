"""
routines_noir.py — NOIR DataLab queries for DESI Legacy Imaging Survey DR10.

Queries the Legacy Survey ``ls_dr10.tractor`` table via the NOIR DataLab
SQL interface.  The ``dl`` library (noirlab-datalab package) must be
installed; if it is not available an ImportError is raised with a clear
installation message.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

try:
    from dl import queryClient as qc
    from dl.helpers.utils import convert
    _DL_AVAILABLE = True
except ImportError:
    _DL_AVAILABLE = False

from astropy import table
import pandas as pd


def _check_dl():
    """Raise ImportError if the noirlab-datalab library is not installed."""
    if not _DL_AVAILABLE:
        raise ImportError(
            "The NOIR DataLab client library is required for Legacy Survey queries.\n"
            "Install it with:  pip install noirlab-datalab")


def query_noir_datalab_stars(RA, DEC, RADIUS=16):
    """Query Legacy Survey DR10 for PSF-type (stellar) sources.

    Parameters
    ----------
    RA : float
        Right ascension in decimal degrees.
    DEC : float
        Declination in decimal degrees.
    RADIUS : float
        Search radius in arcminutes (default: 16).

    Returns
    -------
    astropy.table.Table
        Columns: ra, dec, mag_g, mag_r, mag_i, mag_z,
        snr_g, snr_r, snr_i, snr_z,
        allmask_g, allmask_r, allmask_i, allmask_z.

    Raises
    ------
    ImportError
        If noirlab-datalab is not installed.
    """
    _check_dl()

    query = (
        "SELECT ra, dec, "
        "mag_g, mag_r, mag_i, mag_z, "
        "snr_g, snr_r, snr_i, snr_z, "
        "allmask_g, allmask_r, allmask_i, allmask_z "
        "FROM ls_dr10.tractor AS tractor "
        f"WHERE 't' = q3c_join({RA}, {DEC}, "
        f"tractor.ra, tractor.dec, {RADIUS:.1f}/60.0) "
        "AND tractor.type = 'PSF'"
    )
    result = convert(qc.query(sql=query, timeout=30), 'pandas')
    return table.Table.from_pandas(result)


def query_noir_datalab_extended(RA, DEC, RADIUS=16):
    """Query Legacy Survey DR10 for extended (non-PSF) sources.

    Parameters
    ----------
    RA : float
        Right ascension in decimal degrees.
    DEC : float
        Declination in decimal degrees.
    RADIUS : float
        Search radius in arcminutes (default: 16).

    Returns
    -------
    astropy.table.Table
        Same columns as :func:`query_noir_datalab_stars`.

    Raises
    ------
    ImportError
        If noirlab-datalab is not installed.
    """
    _check_dl()

    query = (
        "SELECT ra, dec, "
        "mag_g, mag_r, mag_i, mag_z, "
        "snr_g, snr_r, snr_i, snr_z, "
        "allmask_g, allmask_r, allmask_i, allmask_z "
        "FROM ls_dr10.tractor AS tractor "
        f"WHERE 't' = q3c_join({RA}, {DEC}, "
        f"tractor.ra, tractor.dec, {RADIUS:.1f}/60.0) "
        "AND NOT tractor.type = 'PSF'"
    )
    result = convert(qc.query(sql=query, timeout=30), 'pandas')
    return table.Table.from_pandas(result)
