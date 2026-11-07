"""Tools for comparing light curve detrending methods on a known transit.

Everything here is plain NumPy, SciPy and pandas so each step can be read
and checked. Time is in days (BTJD). Flux is normalized so the star is 1.0.
"""
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

CADENCE = 2.0 / 1440.0  # 2 minute TESS cadence, in days


# ---------------------------------------------------------------- data ----
def load_data(path):
    """Read the exported CSV. Returns time and the two raw flux columns."""
    d = pd.read_csv(path)
    return d


# --------------------------------------------------------------- shapes ----
def trapezoid_unit(x, t14, t12):
    """Transit shape with depth 1: flat bottom and straight ingress/egress.

    x is time from mid transit, t14 the full duration, t12 the ingress time.
    Returns 1 at the bottom and 0 outside the transit.
    """
    x = np.asarray(x, dtype=float)
    ax = np.abs(x)
    half = t14 / 2.0
    flat_half = half - t12
    s = np.zeros_like(x)
    s[ax <= flat_half] = 1.0
    ramp = (ax > flat_half) & (ax < half)
    s[ramp] = (half - ax[ramp]) / t12
    return s


def time_from_event(t, period, t0, shift=0.0):
    """Time from the nearest event at t0 + (n + shift) * period."""
    return ((t - t0 - shift * period + 0.5 * period) % period) - 0.5 * period


def event_mask(t, period, t0, half_width, shift=0.0):
    """True for points within half_width days of any event."""
    return np.abs(time_from_event(t, period, t0, shift)) < half_width


def inject_transit(t, f, tc, depth, t14, t12):
    """Multiply a transit of known depth into the flux, centered at tc."""
    s = trapezoid_unit(t - tc, t14, t12)
    return f * (1.0 - depth * s)


# ------------------------------------------------------ global detrending ----
def global_median(t, f, mask):
    """No detrending: divide by the median of all unmasked points."""
    return f / np.median(f[~mask])


def split_segments(t, gap_days=0.1):
    """Index ranges of continuous stretches (the orbits), split at gaps."""
    breaks = np.where(np.diff(t) > gap_days)[0] + 1
    starts = np.r_[0, breaks]
    ends = np.r_[breaks, len(t)]
    return list(zip(starts, ends))


def savgol_flatten(t, f, mask, window_days, polyorder=2, gap_days=0.1):
    """Savitzky-Golay trend fitted with masked points filled in.

    Masked points (transits, eclipses) are replaced by a straight line
    between the nearest unmasked points before smoothing, so the smoother
    does not follow the dips. Each orbit is smoothed on its own.
    """
    trend = np.empty_like(f)
    for a, b in split_segments(t, gap_days):
        ts, fs, ms = t[a:b], f[a:b], mask[a:b]
        if (~ms).sum() < 10:
            trend[a:b] = np.median(fs)
            continue
        filled = fs.copy()
        filled[ms] = np.interp(ts[ms], ts[~ms], fs[~ms])
        n = int(round(window_days / CADENCE))
        n = min(n, len(ts))          # cannot exceed the segment
        if n % 2 == 0:
            n -= 1                   # window length must be odd
        n = max(n, 5)
        trend[a:b] = savgol_filter(filled, n, polyorder, mode="interp")
    return f / trend


# ------------------------------------------------------- local detrending ----
def window_info(t, tc, half_window, mask_half, central_half, extra_mask=None):
    """Select the points near one transit and say whether the window is usable."""
    x = t - tc
    in_win = np.abs(x) < half_window
    central = np.abs(x) < central_half
    base = in_win & (np.abs(x) >= mask_half)
    if extra_mask is not None:
        base &= ~extra_mask
    n_exp = 2.0 * central_half / CADENCE
    ok = (
        (in_win & central).sum() >= 0.9 * n_exp
        and (base & (x < 0)).sum() >= 40
        and (base & (x > 0)).sum() >= 40
    )
    return in_win, central, base, ok


def local_poly_detrend(t, f, tc, half_window, mask_half, central_half, degree, extra_mask=None):
    """Fit a polynomial to the baseline around one transit and divide it out.

    Returns the ratio flux / fit for the points in the window (and the window
    mask) or None if the window is not usable.
    """
    in_win, central, base, ok = window_info(t, tc, half_window, mask_half, central_half, extra_mask)
    if not ok:
        return None
    x = (t - tc) / half_window
    coef = np.polynomial.polynomial.polyfit(x[base], f[base], degree)
    trend = np.polynomial.polynomial.polyval(x[in_win], coef)
    return in_win, f[in_win] / trend


# ------------------------------------------------------ depth estimators ----
def depth_central(x, r, central_half):
    """Primary estimator: 1 minus the mean flux in the central part of the transit."""
    sel = np.abs(x) < central_half
    return 1.0 - np.mean(r[sel])


def depth_shape(x, r, t14, t12):
    """Secondary estimator: best scale of the trapezoid shape (least squares)."""
    s = trapezoid_unit(x, t14, t12)
    return np.sum((1.0 - r) * s) / np.sum(s * s)


# ------------------------------------------------------- method registry ----
def make_methods():
    """The list of detrending methods compared in this project."""
    methods = []
    methods.append(dict(name="SAP none", flux="sap", kind="none"))
    for deg in (0, 1, 2, 3):
        methods.append(dict(name=f"SAP local poly {deg}", flux="sap", kind="local", degree=deg))
    for w in (0.5, 1.0, 2.0):
        methods.append(dict(name=f"SAP Savgol {w} d", flux="sap", kind="savgol", window=w))
    methods.append(dict(name="PDCSAP none", flux="pdcsap", kind="none"))
    methods.append(dict(name="PDCSAP local poly 1", flux="pdcsap", kind="local", degree=1))
    return methods


def detrend_global(method, t, f, all_mask):
    """Run a global method (no detrending, or Savitzky-Golay) on the whole curve."""
    if method["kind"] == "none":
        return global_median(t, f, all_mask)
    if method["kind"] == "savgol":
        return savgol_flatten(t, f, all_mask, method["window"])
    raise ValueError("not a global method")


def measure_one(method, t, f, tc, settings, all_mask, r_all=None):
    """Detrend one transit with one method and measure its depth.

    all_mask marks every transit-like event (real, secondary eclipse, injected).
    For global methods, r_all can hold the already detrended full curve so it
    is computed once. Returns a dict, or None when the window is not usable.
    """
    W = settings["half_window"]
    mh = settings["mask_half"]
    ch = settings["central_half"]
    t14, t12 = settings["t14"], settings["t12"]
    x_all = t - tc

    # same usability rule for every method so all compare the same transits
    in_win, central, base, ok = window_info(t, tc, W, mh, ch, extra_mask=all_mask)
    if not ok:
        return None

    if method["kind"] == "local":
        out = local_poly_detrend(t, f, tc, W, mh, ch, method["degree"], extra_mask=all_mask)
        if out is None:
            return None
        win, r = out
        x = x_all[win]
    else:
        if r_all is None:
            r_all = detrend_global(method, t, f, all_mask)
        r = r_all[in_win]
        x = x_all[in_win]

    return dict(
        depth_central=depth_central(x, r, ch),
        depth_shape=depth_shape(x, r, t14, t12),
    )


# --------------------------------------------------------- ephemeris fit ----
def trapezoid_model(params, t):
    """Periodic trapezoid transit: params = period, t0, depth, t14, t12, baseline."""
    period, t0, depth, t14, t12, base = params
    x = time_from_event(t, period, t0)
    return base * (1.0 - depth * trapezoid_unit(x, t14, t12))


def fit_ephemeris(t, y, start):
    """Least squares fit of the periodic trapezoid. Used for timing and shape only."""
    from scipy.optimize import least_squares

    lo = [start[0] - 0.01, start[1] - 0.1, 0.001, 0.05, 0.003, 0.98]
    hi = [start[0] + 0.01, start[1] + 0.1, 0.03, 0.14, 0.05, 1.02]
    fit = least_squares(
        lambda p: y - trapezoid_model(p, t), start, bounds=(lo, hi),
        x_scale=[1e-4, 1e-3, 1e-3, 1e-2, 1e-2, 1e-3],
    )
    resid = y - trapezoid_model(fit.x, t)
    cov = np.linalg.inv(fit.jac.T @ fit.jac) * resid.var()
    return fit.x, np.sqrt(np.diag(cov)), resid.std()


# ------------------------------------------------------- phase curve ----
def phase_curve_design(t, period, t0):
    """Design matrix of a constant plus two harmonics of the orbital period."""
    phi = ((t - t0) / period) % 1.0
    return np.column_stack([
        np.ones_like(t),
        np.cos(2 * np.pi * phi), np.sin(2 * np.pi * phi),
        np.cos(4 * np.pi * phi), np.sin(4 * np.pi * phi),
    ])


def fit_phase_curve(t, y, mask, period, t0):
    """Smooth brightness change over the orbit, fitted outside transit and eclipse.

    Returns the model at every time (normalized to a mean of 1) and the coefficients.
    """
    A = phase_curve_design(t, period, t0)
    coef, *_ = np.linalg.lstsq(A[~mask], y[~mask], rcond=None)
    model = A @ coef
    return model / np.median(model), coef
