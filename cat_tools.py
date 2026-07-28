"""
cat_tools.py — Photometric catalogue metadata, colour transformations,
               positional cross-matching, and Vizier / NOIR DataLab query wrappers.

Cross-matching uses astropy.coordinates.SkyCoord.match_to_catalog_sky()
(C-compiled KD-tree, RA wrap-around safe, ~5-10× faster than the previous
manual KD-tree for 10 000-entry catalogues).

Colour transformation functions follow PEP 8 snake_case naming.  The
Observatory / filter-system naming in function names (e.g. ``ps1``, ``sdss``,
``bessel``) is kept in lower-case to distinguish from all-caps constants.

NOIR DataLab queries (formerly in routines_noir.py) are included here; the
``dl`` library is an optional dependency — an ImportError with a clear
installation message is raised only when a query function is actually called.
"""

from __future__ import annotations

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"


import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy import coordinates as coord, table, units as u
from astropy.io import ascii
from astroquery.vizier import Vizier

from utils import bcolors


# ---------------------------------------------------------------------------
# Catalogue zero-points and properties
# ---------------------------------------------------------------------------

catalog_zp: dict[str, float] = {
    'WISE_W1':   20.500 + 2.699,
    'WISE_W2':   19.500 + 3.339,
    'WISE_W3':   18.000 + 5.174,
    'WISE_W4':   13.000 + 6.620,
    'GALEX_FUV': 18.82,
    'GALEX_NUV': 20.08,
}

catalog_prop: dict[str, Any] = {}

catalog_prop['2MASS'] = {
    'KEYWORDS':   ['RAJ2000', 'DEJ2000', 'Jmag', 'e_Jmag',
                   'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag'],
    'CATID':      'II/246',
    'CATID_OUT':  'II/246/out',
    'FILTER':     ['J', 'H', 'K'],
    'SIGMA_HIGH': 0,
    'SIGMA_LOW':  0.3,
}

catalog_prop['SDSS'] = {
    'KEYWORDS':   ['RA_ICRS', 'DE_ICRS', 'class', 'mode',
                   'umag', 'e_umag', 'gmag', 'e_gmag', 'rmag', 'e_rmag',
                   'imag', 'e_imag', 'zmag', 'e_zmag'],
    'CATID':      'V/147',
    'CATID_OUT':  'V/147/sdss12',
    'FILTER':     ['u', 'g', 'r', 'i', 'z'],
    'SIGMA_HIGH': 0,
    'SIGMA_LOW':  0.2,
}

catalog_prop['PanSTARRS'] = {
    'KEYWORDS':   ['RAJ2000', 'DEJ2000',
                   'o_gmag', 'gPSFf', 'gmag', 'e_gmag',
                   'o_rmag', 'rPSFf', 'rmag', 'e_rmag',
                   'o_imag', 'iPSFf', 'imag', 'e_imag', 'iKmag',
                   'o_zmag', 'zPSFf', 'zmag', 'e_zmag',
                   'o_ymag', 'yPSFf', 'ymag', 'e_ymag'],
    'CATID':      'II/349',
    'CATID_OUT':  'II/349/ps1',
    'FILTER':     ['g', 'r', 'i', 'z', 'y'],
    'SIGMA_HIGH': 0,
    'SIGMA_LOW':  0.2,
}

catalog_prop['GAIA'] = {
    'KEYWORDS':   ['RAJ2000', 'DEJ2000', 'Plx', 'e_Plx',
                   'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE'],
    'CATID':      'I/345',
    'CATID_OUT':  'I/345/gaia2',
}

catalog_prop['UKIDSS'] = {
    'KEYWORDS':   ['RAJ2000', 'DEJ2000', 'class',
                   'Ymag', 'e_Ymag', 'Jmag1', 'e_Jmag1', 'Jmag2', 'e_Jmag2',
                   'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag'],
    'CATID':      'II/319',
    'CATID_OUT':  'II/319/las9',
    'FILTER':     ['Y', 'J', 'H', 'K'],
    'SIGMA_HIGH': 0,
    'SIGMA_LOW':  0.2,
}

catalog_prop['WISE'] = {
    'KEYWORDS':   ['RAJ2000', 'DEJ2000', 'W1mag', 'e_W1mag',
                   'W2mag', 'e_W2mag'],
    'CATID':      'II/328',
    'CATID_OUT':  'II/328/allwise',
    'FILTER':     ['W1', 'W2'],
    'SIGMA_HIGH': 0,
    'SIGMA_LOW':  0.2,
}

catalog_prop['DES'] = {
    'KEYWORDS':   ['RAJ2000', 'DEJ2000',
                   'gmag', 'rmag', 'imag', 'zmag', 'Ymag',
                   'e_gmag', 'e_rmag', 'e_imag', 'e_zmag', 'e_Ymag',
                   'S_Gg', 'S_Gr', 'S_Gi', 'S_Gz',
                   'gFlag', 'rFlag', 'iFlag', 'zFlag'],
    'CATID':      'II/357',
    'CATID_OUT':  'II/357/des_dr1',
    'FILTER':     ['g', 'r', 'i', 'z', 'y'],
    'SIGMA_HIGH': 0,
    'SIGMA_LOW':  0.2,
}


# ---------------------------------------------------------------------------
# Colour transformation functions
# ---------------------------------------------------------------------------

def ps1_to_hsc(data: table.Table) -> table.Table:
    """Convert PanSTARRS photometry to HSC system.

    Reference: https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_8_e/colorterms.html

    Parameters
    ----------
    data : astropy.table.Table
        Table with columns ``{g,r,i,z,y}_PS1`` and ``{g,r,i,z,y}_PS1_ERR``.

    Returns
    -------
    astropy.table.Table
        Input table with added ``{g,r,r2,i,i2,z,y}_HSC`` columns.
    """
    gr = data['g_PS1'] - data['r_PS1']
    ri = data['r_PS1'] - data['i_PS1']
    iz = data['i_PS1'] - data['z_PS1']
    zy = data['z_PS1'] - data['y_PS1']
    yz = data['y_PS1'] - data['z_PS1']

    gr_e = np.sqrt(data['g_PS1_ERR']**2 + data['r_PS1_ERR']**2)
    ri_e = np.sqrt(data['r_PS1_ERR']**2 + data['i_PS1_ERR']**2)
    iz_e = np.sqrt(data['i_PS1_ERR']**2 + data['z_PS1_ERR']**2)
    zy_e = np.sqrt(data['z_PS1_ERR']**2 + data['y_PS1_ERR']**2)
    yz_e = np.sqrt(data['y_PS1_ERR']**2 + data['z_PS1_ERR']**2)

    data['g_HSC']   = data['g_PS1'] + 0.00573 + 0.06175 * gr - 0.001125 * gr**2
    data['r_HSC']   = data['r_PS1'] - 0.00014 + 0.00137 * ri - 0.008380 * ri**2
    data['r2_HSC']  = data['r_PS1'] - 0.000032 - 0.002866 * ri - 0.012638 * ri**2
    data['i_HSC']   = data['i_PS1'] + 0.000643 - 0.130078 * iz - 0.006855 * iz**2
    data['i2_HSC']  = data['i_PS1'] + 0.001625 - 0.200406 * iz - 0.013666 * iz**2
    data['z_HSC']   = data['z_PS1'] - 0.005362 - 0.221551 * zy - 0.308279 * zy**2
    data['y_HSC']   = data['y_PS1'] - 0.002055 + 0.209680 * yz + 0.227296 * yz**2

    data['g_HSC_ERR']  = np.sqrt(data['g_PS1_ERR']**2 + (0.06175*gr_e)**2 + (2*0.001125*gr*gr_e)**2)
    data['r_HSC_ERR']  = np.sqrt(data['r_PS1_ERR']**2 + (0.00137*ri_e)**2 + (2*0.008380*ri*ri_e)**2)
    data['r2_HSC_ERR'] = np.sqrt(data['r_PS1_ERR']**2 + (0.002866*ri_e)**2 + (2*0.012638*ri*ri_e)**2)
    data['i_HSC_ERR']  = np.sqrt(data['i_PS1_ERR']**2 + (0.130078*iz_e)**2 + (2*0.006855*iz*iz_e)**2)
    data['i2_HSC_ERR'] = np.sqrt(data['i_PS1_ERR']**2 + (0.200406*iz_e)**2 + (2*0.013666*iz*iz_e)**2)
    data['z_HSC_ERR']  = np.sqrt(data['z_PS1_ERR']**2 + (0.221551*zy_e)**2 + (2*0.308279*zy*zy_e)**2)
    data['y_HSC_ERR']  = np.sqrt(data['y_PS1_ERR']**2 + (0.209680*yz_e)**2 + (2*0.227296*yz*yz_e)**2)

    return data


def ps1_to_sdss(data: table.Table) -> table.Table:
    """Convert PanSTARRS to SDSS photometry.

    Reference: Finkbeiner et al. (2016), ApJ 822, 66.
    Valid for main-sequence stars with 0.4 < g_PS1 - i_PS1 < 2.7.

    Parameters
    ----------
    data : astropy.table.Table
        Table with PS1 griz y columns.

    Returns
    -------
    astropy.table.Table
        Input table filtered to valid colour range, with added SDSS columns.
    """
    coeffs = {
        'u': [ 0.04438, -2.26095, -0.13387,  0.27099],
        'g': [-0.01808, -0.13595,  0.01941, -0.00183],
        'r': [-0.01836, -0.03577,  0.02612, -0.00558],
        'i': [ 0.01170, -0.00400,  0.00066, -0.00058],
        'z': [-0.01062,  0.07529, -0.03592,  0.00890],
        'y': [ 0.08924, -0.20878,  0.10360, -0.02441],
    }

    x     = data['g_PS1'] - data['i_PS1']
    x_err = np.sqrt(data['g_PS1_ERR']**2 + data['i_PS1_ERR']**2)
    mask  = (x >= 0.4) & (x <= 2.7)
    data  = data[mask]
    x     = x[mask]
    x_err = x_err[mask]

    for filt in ('u', 'g', 'r', 'i', 'z'):
        c     = coeffs[filt]
        poly  = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
        dpoly = np.sqrt((c[1]*x_err)**2 + (2*c[2]*x*x_err)**2 + (3*c[3]*x**2*x_err)**2)
        base  = data['g_PS1'] if filt == 'u' else data[f'{filt}_PS1']
        base_e= data['g_PS1_ERR'] if filt == 'u' else data[f'{filt}_PS1_ERR']
        data[f'{filt}_SDSS']     = np.array(base - poly)
        data[f'{filt}_SDSS_ERR'] = np.array(np.sqrt(base_e**2 + dpoly**2))

    return data


def ps1_to_ztf(data: table.Table) -> table.Table:
    """Convert PanSTARRS to ZTF photometry.

    Reference: Medford et al. (2020), RNAAS 4, 28.

    Parameters
    ----------
    data : astropy.table.Table
        Table with PS1 gr columns.

    Returns
    -------
    astropy.table.Table
        Input table with added ``g_ZTF``, ``r_ZTF`` columns.
    """
    gr = data['g_PS1'] - data['r_PS1']
    data['g_ZTF']     = data['g_PS1'] + 0.055 * gr - 0.012
    data['r_ZTF']     = data['r_PS1'] - 0.087 * gr - 0.0035
    data['g_ZTF_ERR'] = np.sqrt(data['g_PS1_ERR']**2 + (0.055*data['g_PS1_ERR'])**2 + (0.055*data['r_PS1_ERR'])**2)
    data['r_ZTF_ERR'] = np.sqrt(data['r_PS1_ERR']**2 + (0.087*data['g_PS1_ERR'])**2 + (0.087*data['r_PS1_ERR'])**2)
    return data


# PS1 grizy (AB) -> 2MASS J (Vega) synthetic transformation.
# Derived by synthetic photometry of the Pickles stellar atlas through the PS1
# and 2MASS passbands (SVO) with Fitzpatrick (1999) reddening; see the notebook
# ps1_to_2mass_J_transformation.ipynb.  Cross-checked against Tonry et al.
# (2012): synthetic A0V y_P1 - J = 0.530 vs their tabulated 0.531.
#     J_2MASS = y_P1 - P(g-i),   P a cubic in (g-i)_P1
PS1_J_COEFFS   = np.array([-0.07038, 0.01407, 0.40652, 0.59905])  # numpy high->low order
PS1_J_GI_RANGE = (-0.86, 1.69)   # applicability: excludes M dwarfs+ (red) and very blue (O/early-B)
PS1_J_SYS      = 0.05            # synthetic transformation scatter (mag)


def ps1_to_2mass_j(data: table.Table) -> table.Table:
    """Predict 2MASS ``J`` (Vega) from PanSTARRS1 ``grizy`` (AB) for stars.

    ``J = y_PS1 - P(g-i)`` with ``P`` the cubic :data:`PS1_J_COEFFS`, derived
    from synthetic photometry of the Pickles atlas through the PS1 and 2MASS
    passbands (see ``ps1_to_2mass_J_transformation.ipynb``).  Valid for roughly
    B--K stars, ``-0.86 < (g-i)_PS1 < 1.69``; stars outside that range
    (M dwarfs and later, and very blue O/early-B stars) are set to NaN so that
    they are dropped downstream.

    Parameters
    ----------
    data : astropy.table.Table
        Table with PS1 ``g_PS1``, ``i_PS1``, ``y_PS1`` columns and their
        ``*_PS1_ERR`` uncertainties.

    Returns
    -------
    astropy.table.Table
        Input table with added ``J_2MASS`` (Vega mag) and ``J_2MASS_ERR``.
    """
    gi = np.asarray(data['g_PS1'] - data['i_PS1'], dtype=float)
    P  = np.polyval(PS1_J_COEFFS, gi)
    dP = np.polyval(np.polyder(PS1_J_COEFFS), gi)          # dP/d(g-i)
    J  = np.asarray(data['y_PS1'], float) - P
    # sigma_J^2 = sigma_y^2 + (dP)^2 (sigma_g^2 + sigma_i^2) + sigma_sys^2
    Jerr = np.sqrt(np.asarray(data['y_PS1_ERR'], float)**2
                   + dP**2 * (np.asarray(data['g_PS1_ERR'], float)**2
                              + np.asarray(data['i_PS1_ERR'], float)**2)
                   + PS1_J_SYS**2)
    outside = (gi < PS1_J_GI_RANGE[0]) | (gi > PS1_J_GI_RANGE[1]) | ~np.isfinite(gi)
    data['J_2MASS']     = np.where(outside, np.nan, J)
    data['J_2MASS_ERR'] = np.where(outside, np.nan, Jerr)
    return data


def stellar_locus_mask(data: table.Table, nsigma: float = 3.0,
                       niter: int = 3) -> np.ndarray:
    """Flag PS1 stars lying off the fiducial ``(g-r, r-i)`` stellar locus.

    A robust cubic is fit to ``(r-i)`` vs ``(g-r)`` and outliers are rejected
    iteratively at *nsigma* times the MAD-based scatter.  This statistically
    removes strongly non-solar-metallicity stars, unresolved binaries, blends,
    QSO/galaxy contaminants and bad photometry, which lie off the tight
    dwarf/giant locus.  Reddening moves stars *along* the locus, so reddened
    stars are not preferentially rejected.

    Note: broadband PS1 colours cannot separate luminosity class, so normal-
    coloured supergiants are *not* removed by this cut — they are rare and get
    down-weighted / clipped in the zeropoint bootstrap.  A Gaia parallax or
    reduced-proper-motion cut would be needed to target them directly.

    Parameters
    ----------
    data : astropy.table.Table
        Table with PS1 ``g_PS1``, ``r_PS1``, ``i_PS1`` columns.
    nsigma : float
        Rejection threshold in robust sigma.
    niter : int
        Number of clipping iterations.

    Returns
    -------
    numpy.ndarray of bool
        True for on-locus stars.
    """
    gr = np.asarray(data['g_PS1'] - data['r_PS1'], float)
    ri = np.asarray(data['r_PS1'] - data['i_PS1'], float)
    finite = np.isfinite(gr) & np.isfinite(ri)
    keep = finite.copy()
    for _ in range(niter):
        if keep.sum() < 8:
            break
        p     = np.polyfit(gr[keep], ri[keep], 3)
        resid = ri - np.polyval(p, gr)
        med   = np.median(resid[keep])
        mad   = np.median(np.abs(resid[keep] - med))
        sig   = 1.4826 * mad if mad > 0 else np.std(resid[keep])
        if sig <= 0:
            break
        keep = finite & (np.abs(resid - med) < nsigma * sig)
    return keep


def sdss_to_des(data: table.Table) -> table.Table:
    """Convert SDSS to DES photometry.

    Reference: Drlica-Wagner et al. (2018), ApJS 235, 33, Appendix A.4.

    Parameters
    ----------
    data : astropy.table.Table
        Table with SDSS griz columns.

    Returns
    -------
    astropy.table.Table
        Input table with added DES griz columns.
    """
    data['g_DES'] = data['g_SDSS'] - 0.104*(data['g_SDSS'] - data['r_SDSS']) + 0.01
    data['r_DES'] = data['r_SDSS'] - 0.102*(data['g_SDSS'] - data['r_SDSS']) + 0.02
    data['i_DES'] = data['i_SDSS'] - 0.256*(data['i_SDSS'] - data['z_SDSS']) + 0.02
    data['z_DES'] = data['z_SDSS'] - 0.086*(data['i_SDSS'] - data['z_SDSS']) + 0.01

    data['g_DES_ERR'] = np.sqrt(data['g_SDSS_ERR']**2 + 0.104**2*(data['g_SDSS_ERR']**2 + data['r_SDSS_ERR']**2))
    data['r_DES_ERR'] = np.sqrt(data['r_SDSS_ERR']**2 + 0.102**2*(data['g_SDSS_ERR']**2 + data['r_SDSS_ERR']**2))
    data['i_DES_ERR'] = np.sqrt(data['i_SDSS_ERR']**2 + 0.256**2*(data['i_SDSS_ERR']**2 + data['z_SDSS_ERR']**2))
    data['z_DES_ERR'] = np.sqrt(data['z_SDSS_ERR']**2 + 0.086**2*(data['i_SDSS_ERR']**2 + data['z_SDSS_ERR']**2))

    return data


def sdss_to_grond(data: table.Table) -> table.Table:
    """Convert SDSS to GROND photometry.

    Reference: https://www.mpe.mpg.de/~jcg/GROND/calibration.html

    Parameters
    ----------
    data : astropy.table.Table
        Table with SDSS griz columns.

    Returns
    -------
    astropy.table.Table
        Input table with added GROND griz columns.
    """
    data['g_GROND'] = data['g_SDSS'] - 0.015*(data['g_SDSS'] - data['r_SDSS']) + 0.006
    data['r_GROND'] = data['r_SDSS'] - 0.012*(data['r_SDSS'] - data['i_SDSS']) + 0.004
    data['i_GROND'] = data['i_SDSS'] - 0.113*(data['r_SDSS'] - data['i_SDSS']) + 0.031
    data['z_GROND'] = data['z_SDSS'] + 0.009*(data['i_SDSS'] - data['z_SDSS']) + 0.003

    data['g_GROND_ERR'] = np.sqrt(data['g_SDSS_ERR']**2 + 0.015**2*(data['g_SDSS_ERR']**2 + data['r_SDSS_ERR']**2))
    data['r_GROND_ERR'] = np.sqrt(data['r_SDSS_ERR']**2 + 0.012**2*(data['r_SDSS_ERR']**2 + data['i_SDSS_ERR']**2))
    data['i_GROND_ERR'] = np.sqrt(data['i_SDSS_ERR']**2 + 0.113**2*(data['r_SDSS_ERR']**2 + data['i_SDSS_ERR']**2))
    data['z_GROND_ERR'] = np.sqrt(data['z_SDSS_ERR']**2 + 0.009**2*(data['i_SDSS_ERR']**2 + data['z_SDSS_ERR']**2))

    return data


def sdss_to_bessel(data: table.Table) -> table.Table:
    """Convert SDSS to Bessel (Johnson-Cousins) photometry.

    Reference: http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php

    Parameters
    ----------
    data : astropy.table.Table
        Table with SDSS griz columns.

    Returns
    -------
    astropy.table.Table
        Input table with added Bessel BVRI columns.
    """
    data['B_BESSEL'] = data['g_SDSS'] + 0.3130*(data['g_SDSS'] - data['r_SDSS']) + 0.2271
    data['B_BESSEL_ERR'] = np.sqrt(0.0107**2 + data['g_SDSS_ERR']**2
                                    + (0.3130*data['g_SDSS_ERR'])**2
                                    + (0.3130*data['r_SDSS_ERR'])**2)

    data['V_BESSEL'] = data['g_SDSS'] - 0.5784*(data['g_SDSS'] - data['r_SDSS']) - 0.0038
    data['V_BESSEL_ERR'] = np.sqrt(0.0054**2 + data['g_SDSS_ERR']**2
                                    + (0.5784*data['g_SDSS_ERR'])**2
                                    + (0.5784*data['r_SDSS_ERR'])**2)

    data['R_BESSEL'] = data['r_SDSS'] - 0.1837*(data['g_SDSS'] - data['r_SDSS']) - 0.0971
    data['R_BESSEL_ERR'] = np.sqrt(0.0106**2 + data['r_SDSS_ERR']**2
                                    + (0.1837*data['g_SDSS_ERR'])**2
                                    + (0.1837*data['r_SDSS_ERR'])**2)

    data['I_BESSEL'] = data['r_SDSS'] - 1.2444*(data['r_SDSS'] - data['i_SDSS']) - 0.3820
    data['I_BESSEL_ERR'] = np.sqrt(0.0078**2 + data['r_SDSS_ERR']**2
                                    + (1.2444*data['r_SDSS_ERR'])**2
                                    + (1.2444*data['i_SDSS_ERR'])**2)

    return data


def des_to_sdss(data: table.Table) -> table.Table:
    """Convert DES to SDSS photometry (inverted transformation).

    Reference: Drlica-Wagner et al. (2018), Appendix A.4, inverted with Mathematica.

    Parameters
    ----------
    data : astropy.table.Table
        Table with DES griz columns.

    Returns
    -------
    astropy.table.Table
        Input table with added SDSS griz columns.
    """
    data['g_SDSS'] =  1.002  * (-0.00894 + 1.102*data['g_DES'] - 0.104*data['r_DES'])
    data['r_SDSS'] = -1.002  * ( 0.01894 - 0.102*data['g_DES'] - 0.896*data['r_DES'])
    data['i_SDSS'] = -1.20482* ( 0.01916 - 1.086*data['i_DES'] + 0.256*data['z_DES'])
    data['z_SDSS'] = -1.20482* ( 0.00916 - 0.086*data['i_DES'] - 0.744*data['z_DES'])

    data['g_SDSS_ERR'] = np.sqrt((1.002  *1.102*data['g_DES_ERR'])**2 + (1.002  *0.104*data['r_DES_ERR'])**2)
    data['r_SDSS_ERR'] = np.sqrt((1.002  *0.102*data['g_DES_ERR'])**2 + (1.002  *0.896*data['r_DES_ERR'])**2)
    data['i_SDSS_ERR'] = np.sqrt((1.20482*1.086*data['i_DES_ERR'])**2 + (1.20482*0.256*data['z_DES_ERR'])**2)
    data['z_SDSS_ERR'] = np.sqrt((1.20482*0.086*data['i_DES_ERR'])**2 + (1.20482*0.744*data['z_DES_ERR'])**2)

    return data


# ---------------------------------------------------------------------------
# Vizier query wrapper
# ---------------------------------------------------------------------------

def retrieve_photcat(
    object_prop: dict[str, Any],
    photcat: str,
    catprop: dict[str, Any],
    filename: str | Path | None = None,
    row_limit: int = -1,
    radius: u.Quantity = 10.0 * u.arcmin,
    output_dir: str | Path = 'photcat/',
) -> None:
    """Query Vizier for a photometric catalogue and write matched sources to file.

    Parameters
    ----------
    object_prop : dict
        Must contain ``'RA'`` (decimal deg), ``'DEC'`` (decimal deg),
        ``'FILTER'`` (list of filter strings).
    photcat : str
        Catalogue key, e.g. ``'SDSS'``, ``'PanSTARRS'``, ``'2MASS'``.
    catprop : dict
        Catalogue property dictionary — normally ``cat_tools.catalog_prop``.
    filename : str | Path | None
        Output file path.
    row_limit : int
        Vizier row limit (-1 = unlimited, default).
    radius : astropy.units.Quantity
        Search radius (default: 10 arcmin).
    output_dir : str | Path
        Output directory, created if missing.
    """
    v      = Vizier(columns=['all'], row_limit=row_limit)
    result = v.query_region(
        coord.SkyCoord(object_prop['RA'], object_prop['DEC'], unit=u.deg),
        radius=radius,
        catalog=catprop[photcat]['CATID'])
    result = result[catprop[photcat]['CATID_OUT']]
    result = result[catprop[photcat]['KEYWORDS']]

    for key in [k for k in result.colnames if 'mag' in k]:
        result[key].format = '.4f'

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if photcat == 'SDSS':
        result = result[(result['class'] == 6) & (result['mode'] == 1)]

    if photcat == 'PanSTARRS':
        result = result[
            (result['o_gmag'] >= 0.85) & (result['o_rmag'] >= 0.85) &
            (result['o_imag'] >= 0.85) & (result['o_zmag'] >= 0.85) &
            (result['o_ymag'] >= 0.85) &
            (np.round(result['gPSFf'], 0) >= 1) &
            (np.round(result['rPSFf'], 0) >= 1) &
            (np.round(result['iPSFf'], 0) >= 1) &
            (np.round(result['zPSFf'], 0) >= 1) &
            (np.round(result['yPSFf'], 0) >= 1) &
            (np.abs(result['imag'] - result['iKmag']) < 0.05) &
            (result['imag'] > 14) & (result['imag'] < 21)]

    for filt in object_prop['FILTER']:
        if filt not in catprop[photcat]['FILTER']:
            print(bcolors.FAIL
                  + f'Filter {filt} not in catalogue {photcat}.'
                  + bcolors.ENDC)
            sys.exit(1)

        ra_key  = catprop[photcat]['KEYWORDS'][0]
        dec_key = catprop[photcat]['KEYWORDS'][1]
        mask    = ((result[f'e_{filt}mag'] > catprop[photcat]['SIGMA_HIGH'])
                   & (result[f'e_{filt}mag'] < catprop[photcat]['SIGMA_LOW']))
        ascii.write(result[[ra_key, dec_key, f'{filt}mag', f'e_{filt}mag']][mask],
                    str(filename), overwrite=True, format='no_header')


# ---------------------------------------------------------------------------
# Sky cross-matching
# ---------------------------------------------------------------------------

def crossmatch_sky(
    ra1: np.ndarray,
    dec1: np.ndarray,
    ra2: np.ndarray,
    dec2: np.ndarray,
    max_sep_arcsec: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cross-match two sets of sky coordinates using astropy SkyCoord.

    Finds for each point in catalog 1 the nearest neighbour in catalog 2.
    Handles RA wrap-around correctly.

    Parameters
    ----------
    ra1, dec1 : np.ndarray
        Coordinates of the first catalog (decimal degrees).
    ra2, dec2 : np.ndarray
        Coordinates of the second catalog (decimal degrees).
    max_sep_arcsec : float
        Maximum separation to accept as a match (arcseconds).

    Returns
    -------
    idx : np.ndarray[int]
        Index into catalog 2 of the nearest neighbour for each catalog-1 point.
    sep_arcsec : np.ndarray[float]
        Angular separation in arcseconds for each match.
    mask : np.ndarray[bool]
        True where separation < *max_sep_arcsec*.
    """
    cat1 = coord.SkyCoord(ra=np.asarray(ra1)*u.deg, dec=np.asarray(dec1)*u.deg)
    cat2 = coord.SkyCoord(ra=np.asarray(ra2)*u.deg, dec=np.asarray(dec2)*u.deg)
    idx, d2d, _ = cat1.match_to_catalog_sky(cat2)
    sep         = d2d.arcsec
    return idx, sep, sep < max_sep_arcsec


def crossmatch_catalogs(
    catalog1: np.ndarray,
    catalog2: np.ndarray,
    tolerance_arcsec: float,
) -> np.ndarray:
    """Cross-match two float arrays by RA/Dec (first two columns each).

    Parameters
    ----------
    catalog1 : np.ndarray, shape (N1, M1)
        Columns 0 and 1 are RA, Dec in decimal degrees.
    catalog2 : np.ndarray, shape (N2, M2)
        Columns 0 and 1 are RA, Dec in decimal degrees.
    tolerance_arcsec : float
        Maximum matching radius in arcseconds.

    Returns
    -------
    np.ndarray, shape (N_matched, M1 + M2 + 1)
        Matched rows: ``[catalog2_cols | catalog1_cols | sep_deg]``.
        Separation is in **degrees** (multiply by 3600 for arcsec).
    """
    d1 = np.asarray(catalog1, dtype=float)
    d2 = np.asarray(catalog2, dtype=float)

    if not len(d1) or not len(d2):
        empty_w = d1.shape[1] + d2.shape[1] + 1
        return np.empty((0, empty_w))

    idx, sep, mask = crossmatch_sky(
        d1[:, 0], d1[:, 1],
        d2[:, 0], d2[:, 1],
        tolerance_arcsec)

    matched_d2  = d2[idx[mask]]
    matched_d1  = d1[mask]
    sep_deg     = sep[mask] / 3600.0

    return np.hstack([matched_d2, matched_d1, sep_deg.reshape(-1, 1)])


def crossmatch_angular(
    X1: np.ndarray,
    X2: np.ndarray,
    max_distance: float = np.inf,
) -> tuple[np.ndarray, np.ndarray]:
    """Angular cross-match (legacy wrapper around crossmatch_sky).

    Parameters
    ----------
    X1, X2 : np.ndarray, shape (N, 2)
        RA/Dec arrays in degrees.
    max_distance : float
        Maximum separation in degrees.

    Returns
    -------
    dist : np.ndarray
        Angular separations in degrees (inf for unmatched).
    ind : np.ndarray
        Indices into X2.
    """
    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)

    cat1 = coord.SkyCoord(ra=X1[:, 0]*u.deg, dec=X1[:, 1]*u.deg)
    cat2 = coord.SkyCoord(ra=X2[:, 0]*u.deg, dec=X2[:, 1]*u.deg)
    idx, d2d, _ = cat1.match_to_catalog_sky(cat2)

    dist = d2d.deg
    if not np.isinf(max_distance):
        dist[dist > max_distance] = np.inf

    return dist, idx


# Keep old name for any existing calling code
def wrapper_crossmatch(
    FILE1: np.ndarray,
    FILE2: np.ndarray,
    RADIUS: float,
) -> np.ndarray:
    """Backward-compatible alias for :func:`crossmatch_catalogs`."""
    return crossmatch_catalogs(FILE1, FILE2, RADIUS)


# ---------------------------------------------------------------------------
# NOIR DataLab queries (formerly in routines_noir.py)
# ---------------------------------------------------------------------------

def _require_dl() -> tuple[Any, Any]:
    """Import noirlab-datalab or raise a descriptive ImportError."""
    try:
        from dl import queryClient as qc
        from dl.helpers.utils import convert
        return qc, convert
    except ImportError as exc:
        raise ImportError(
            "The NOIR DataLab library is required for Legacy Survey queries.\n"
            "Install with:  pip install noirlab-datalab"
        ) from exc


def query_noir_datalab_stars(
    ra: float,
    dec: float,
    radius: float = 16.0,
) -> table.Table:
    """Query Legacy Survey DR10 for PSF-type (stellar) sources.

    Parameters
    ----------
    ra : float
        Right ascension in decimal degrees.
    dec : float
        Declination in decimal degrees.
    radius : float
        Search radius in arcminutes (default: 16).

    Returns
    -------
    astropy.table.Table
        Columns: ra, dec, mag_{g,r,i,z}, snr_{g,r,i,z}, allmask_{g,r,i,z}.

    Raises
    ------
    ImportError
        If noirlab-datalab is not installed.
    """
    qc, convert = _require_dl()

    sql = (
        "SELECT ra, dec, "
        "mag_g, mag_r, mag_i, mag_z, "
        "snr_g, snr_r, snr_i, snr_z, "
        "allmask_g, allmask_r, allmask_i, allmask_z "
        "FROM ls_dr10.tractor AS t "
        f"WHERE 't' = q3c_join({ra}, {dec}, t.ra, t.dec, {radius:.1f}/60.0) "
        "AND t.type = 'PSF'"
    )
    return table.Table.from_pandas(convert(qc.query(sql=sql, timeout=30), 'pandas'))


def query_noir_datalab_extended(
    ra: float,
    dec: float,
    radius: float = 16.0,
) -> table.Table:
    """Query Legacy Survey DR10 for extended (non-PSF) sources.

    Parameters
    ----------
    ra : float
        Right ascension in decimal degrees.
    dec : float
        Declination in decimal degrees.
    radius : float
        Search radius in arcminutes (default: 16).

    Returns
    -------
    astropy.table.Table
        Same columns as :func:`query_noir_datalab_stars`.
    """
    qc, convert = _require_dl()

    sql = (
        "SELECT ra, dec, "
        "mag_g, mag_r, mag_i, mag_z, "
        "snr_g, snr_r, snr_i, snr_z, "
        "allmask_g, allmask_r, allmask_i, allmask_z "
        "FROM ls_dr10.tractor AS t "
        f"WHERE 't' = q3c_join({ra}, {dec}, t.ra, t.dec, {radius:.1f}/60.0) "
        "AND NOT t.type = 'PSF'"
    )
    return table.Table.from_pandas(convert(qc.query(sql=sql, timeout=30), 'pandas'))
