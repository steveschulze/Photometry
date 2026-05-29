"""
utils.py — Shared utilities for the ENGRAVE photometry pipeline.

Consolidates terminal colour codes (formerly misc.py) and the bootstrap
zeropoint statistics function (formerly stat_tools.statNclip).
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"


import numpy as np


# ---------------------------------------------------------------------------
# Terminal colour codes
# ---------------------------------------------------------------------------

class bcolors:
    """ANSI escape codes for coloured terminal output."""
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'


# ---------------------------------------------------------------------------
# Statistical utilities
# ---------------------------------------------------------------------------

def bootstrap_zp_stats(
    measurements: np.ndarray,
    niter: int = 1000,
) -> np.ndarray:
    """Robust bootstrap statistics for photometric zeropoint arrays.

    Performs IQR-based outlier rejection (1.5 × IQR fence), then draws
    *niter* bootstrap samples with per-measurement Gaussian scatter to
    estimate the median zeropoint and its asymmetric 1σ uncertainties.

    The entire bootstrap ensemble is drawn in two vectorised NumPy calls,
    giving ~100–500× speedup over a Python loop for niter = 1 000–30 000.

    Parameters
    ----------
    measurements : np.ndarray, shape (N, 2)
        Array of zeropoint measurements.
        ``measurements[:, 0]`` are ZP values; ``measurements[:, 1]`` are
        their 1σ uncertainties.
    niter : int
        Number of bootstrap iterations (default: 1 000).

    Returns
    -------
    np.ndarray, shape (4,)
        ``[zp_median, zp_err_plus, zp_err_minus, n_stars]`` where
        *zp_err_plus* and *zp_err_minus* are the distances from the median
        to the 84th and 16th bootstrap percentiles.
    """
    A = np.asarray(measurements, dtype=float)
    if A.ndim != 2 or A.shape[1] < 2:
        raise ValueError("measurements must have shape (N, 2)")

    zp_vals = A[:, 0]
    zp_errs = A[:, 1]

    # IQR-based outlier rejection
    p25, p75 = np.percentile(zp_vals, [25, 75])
    iqr      = p75 - p25
    mask     = (zp_vals > p25 - 1.5 * iqr) & (zp_vals < p75 + 1.5 * iqr)
    vals_ok  = zp_vals[mask]
    errs_ok  = zp_errs[mask]

    n = len(vals_ok)
    if n == 0:
        return np.array([np.nan, np.nan, np.nan, 0.0])

    # Monte Carlo Gaussian scatter: shape (niter, n)
    random_mc = np.random.normal(vals_ok, errs_ok, size=(niter, n))

    if n > 10:
        # Bootstrap resample with replacement inside each MC draw
        boot_idx     = np.random.randint(0, n, size=(niter, n))
        rows         = np.arange(niter)[:, None]
        random_array = np.median(random_mc[rows, boot_idx], axis=1)
    else:
        random_array = np.median(random_mc, axis=1)

    zp_med = float(np.percentile(random_array, 50))
    zp_sup = float(np.percentile(random_array, 50 + 34.1))
    zp_inf = float(np.percentile(random_array, 50 - 34.1))

    return np.array([zp_med, zp_sup - zp_med, zp_med - zp_inf, float(n)])
