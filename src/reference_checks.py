"""Two consistency checks that use published numbers. They are checks, not results.

1. What bias would the orbital brightness change alone cause for a flat or curved
   baseline, at the real transit position and at the injection positions?
2. What central depth would published system parameters predict?

Usage: python src/reference_checks.py <csv> <results folder>
"""
import json, sys
import numpy as np, pandas as pd
import transit_detrend as td

csv_path, outdir = sys.argv[1], sys.argv[2]
S = json.load(open(f"{outdir}/analysis_settings.json"))
P, T0, T14, T12 = S["period"], S["t0"], S["t14"], S["t12"]
SET = dict(half_window=S["half_window"], mask_half=S["mask_half"], central_half=S["central_half"])

d = pd.read_csv(csv_path)
t = d.time_btjd.values
y = d.pdcsap_flux.values / np.median(d.pdcsap_flux.values)
mask = td.event_mask(t, P, T0, SET["mask_half"]) | td.event_mask(t, P, T0, SET["mask_half"], shift=0.5)
model, coef = td.fit_phase_curve(t, y, mask, P, T0)

# ---- check 1: predicted bias from the fitted orbital brightness change alone
tt = np.arange(1360.0, 1366.0, td.CADENCE)
A = td.phase_curve_design(tt, P, T0)
clean = A @ coef
clean = clean / np.median(clean)
D = 0.0105
pred = {}
for label, shift in [("transit position (phase 0)", 0.0), ("injection position (phase 0.25)", 0.25), ("injection position (phase 0.75)", 0.75)]:
    n_first, n_last = int((1361 - T0) / P), int((1365 - T0) / P)
    out = {0: [], 2: []}
    for n in range(n_first, n_last):
        tc = T0 + (n + shift) * P
        f = td.inject_transit(tt, clean, tc, D, T14, T12)
        m_all = (td.event_mask(tt, P, T0, SET["mask_half"]) | td.event_mask(tt, P, T0, SET["mask_half"], shift=0.5)
                 | (np.abs(tt - tc) < SET["mask_half"]))
        for deg in out:
            r = td.local_poly_detrend(tt, f, tc, SET["half_window"], SET["mask_half"], SET["central_half"], deg, extra_mask=m_all)
            if r is None:
                continue
            win, ratio = r
            out[deg].append(td.depth_central(tt[win] - tc, ratio, SET["central_half"]) - D)
    pred[label] = {f"degree_{k}_bias_ppm": float(np.mean(v) * 1e6) for k, v in out.items()}

# ---- check 2: central depth predicted from published parameters
# Shporer et al. 2019, AJ, Table 2 (radius ratio, impact parameter, a/Rs) and quadratic limb darkening
# from Claret 2017 as used there. Small planet on a grid, circular orbit.
k, b, a = 0.09716, 0.318, 3.562
u1, u2 = 0.2192, 0.3127
Ppub = 0.9414526
cosi = b / a
norm = np.pi * (1 - u1 / 3 - u2 / 6)
g = np.linspace(-1, 1, 301)
X, Y = np.meshgrid(g, g)
inside = X ** 2 + Y ** 2 <= 1.0
px, py, area = k * X[inside], k * Y[inside], (k * (g[1] - g[0])) ** 2

def depth_at(x_days):
    ph = 2 * np.pi * x_days / Ppub
    sep = a * np.sqrt(np.sin(ph) ** 2 + (cosi * np.cos(ph)) ** 2)
    rr2 = (px + sep) ** 2 + py ** 2
    on = rr2 < 1.0
    mu = np.sqrt(1 - rr2[on])
    inten = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
    return np.sum(inten) * area / norm

xs = np.linspace(-SET["central_half"], SET["central_half"], 201)
expected = float(np.mean([depth_at(v) for v in xs]) * 1e6)

published = dict(B1_cos1_ppm=-174.4, B2_cos2_ppm=-190.5, A1_sin1_ppm=21.0,
                 transit_depth_rp_over_rs_squared_ppm=9439, source="Shporer et al. 2019, AJ, Table 2")
ours = dict(zip(["cos1_ppm", "sin1_ppm", "cos2_ppm"], [float(coef[1] * 1e6), float(coef[2] * 1e6), float(coef[3] * 1e6)]))
res = dict(predicted_bias_from_orbital_brightness_change=pred,
           expected_central_depth_from_published_parameters_ppm=expected,
           radius_ratio_squared_ppm=k ** 2 * 1e6,
           our_fitted_harmonics=ours, published_harmonics=published)
json.dump(res, open(f"{outdir}/reference_checks.json", "w"), indent=2)
print(json.dumps(res, indent=2))
