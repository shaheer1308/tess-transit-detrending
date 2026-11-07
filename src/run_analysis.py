"""Full analysis: real WASP-18 sector 2 data plus injection and recovery tests.

Run:  python run_analysis.py <path to wasp18_sector2.csv> <output folder>
"""
import json, sys, os
import numpy as np, pandas as pd
import transit_detrend as td

csv_path = sys.argv[1]
outdir = sys.argv[2]
os.makedirs(outdir, exist_ok=True)

# ---------------------------------------------------------------- data ----
d = td.load_data(csv_path)
t = d.time_btjd.values
flux = {
    "sap": d.sap_flux.values / np.median(d.sap_flux.values),
    "pdcsap": d.pdcsap_flux.values / np.median(d.pdcsap_flux.values),
}

# ------------------------------------------------------------ ephemeris ----
# Step 1: Box Least Squares for a first guess. Step 2: trapezoid fit to PDCSAP.
from astropy.timeseries import BoxLeastSquares
periods = np.arange(0.5, 3.0, 0.00005)
bls = BoxLeastSquares(t, flux["pdcsap"], dy=d.pdcsap_flux_err.values / np.median(d.pdcsap_flux.values))
res = bls.power(periods, np.array([0.05, 0.07, 0.08, 0.09, 0.10, 0.11]), objective="likelihood")
ib = int(np.argmax(res.power))
bls_P, bls_T0 = float(res.period[ib]), float(res.transit_time[ib])
params, sig, rms = td.fit_ephemeris(t, flux["pdcsap"], [bls_P, bls_T0, 0.01, 0.09, 0.015, 1.0])
P, T0, depth_fit, T14, T12, base_fit = [float(v) for v in params]

SET = dict(
    half_window=0.17,          # days each side of mid transit
    mask_half=0.75 * T14,      # points inside this are never used for the baseline fit
    central_half=0.02,         # depth is measured inside this half width
    t14=T14, t12=T12,
)
settings_out = dict(
    bls_period=bls_P, bls_t0=bls_T0,
    period=P, t0=T0, t14=T14, t12=T12, trapezoid_depth_pdcsap=depth_fit,
    one_sigma=dict(zip(["period", "t0", "depth", "t14", "t12", "baseline"], [float(v) for v in sig])),
    rms_residual_ppm=float(rms * 1e6),
    **{k: float(v) for k, v in SET.items() if k not in ("t14", "t12")},
    n_points=int(len(t)),
)

transit_mask = td.event_mask(t, P, T0, SET["mask_half"])
eclipse_mask = td.event_mask(t, P, T0, SET["mask_half"], shift=0.5)
real_mask = transit_mask | eclipse_mask

# Smooth orbital brightness change (reflection and ellipsoidal light), fitted to PDCSAP
# outside transit and eclipse. Used only for a check, never for the main results.
phase_model, phase_coef = td.fit_phase_curve(t, flux["pdcsap"], real_mask, P, T0)
flux_variants = {
    "as observed": flux,
    "orbital brightness change removed": {k: v / phase_model for k, v in flux.items()},
}
settings_out["phase_curve_coefficients_ppm"] = dict(
    zip(["constant_minus_1", "cos_1", "sin_1", "cos_2", "sin_2"],
        [float((phase_coef[0] - 1) * 1e6)] + [float(c * 1e6) for c in phase_coef[1:]]))

segments = td.split_segments(t)
seg_start = np.empty_like(t)
for a, b in segments:
    seg_start[a:b] = t[a]
since_start = t - seg_start

methods = td.make_methods()

# ------------------------------------------------------- helper pieces ----
n_lo, n_hi = int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 1
order = [m["name"] for m in methods]


def real_depths(fl, settings=SET):
    """Measure every method on every usable real transit. Returns a long table."""
    cache = {m["name"]: td.detrend_global(m, t, fl[m["flux"]], real_mask)
             for m in methods if m["kind"] in ("none", "savgol")}
    rows, unusable = [], []
    for n in range(n_lo, n_hi + 1):
        tc = T0 + n * P
        if tc < t.min() or tc > t.max():
            continue
        out = {m["name"]: td.measure_one(m, t, fl[m["flux"]], tc, settings, real_mask,
                                          r_all=cache.get(m["name"])) for m in methods}
        if any(v is None for v in out.values()):
            unusable.append((n, tc))
            continue
        i0 = min(np.searchsorted(t, tc), len(t) - 1)
        for name, r in out.items():
            rows.append(dict(n=n, tc=tc, method=name, depth_central=r["depth_central"],
                             depth_shape=r["depth_shape"], since_orbit_start=float(since_start[i0])))
    return pd.DataFrame(rows), unusable


def summarize_real(df, col):
    g = df.groupby("method")[col]
    s = pd.DataFrame({"n_transits": g.count(), "mean_depth": g.mean(),
                      "std_between_transits": g.std(), "std_error": g.std() / np.sqrt(g.count())})
    return s.loc[order]


def paired_vs_sap_none(df):
    piv = df.pivot(index="n", columns="method", values="depth_central")[order]
    diff = piv.sub(piv["SAP none"], axis=0)
    return pd.DataFrame({"mean_diff_vs_SAP_none": diff.mean(),
                         "std_error_of_diff": diff.std() / np.sqrt(len(piv))})


def injections(fl, depths, shifts):
    rows = []
    for depth in depths:
        for shift in shifts:
            for n in range(n_lo, n_hi + 1):
                tc = T0 + (n + shift) * P
                if tc < t.min() or tc > t.max():
                    continue
                all_mask = real_mask | (np.abs(t - tc) < SET["mask_half"])
                results = {}
                for fluxname in ("sap", "pdcsap"):
                    f_inj = td.inject_transit(t, fl[fluxname], tc, depth, T14, T12)
                    for m in methods:
                        if m["flux"] != fluxname:
                            continue
                        r_all = td.detrend_global(m, t, f_inj, all_mask) if m["kind"] in ("none", "savgol") else None
                        results[m["name"]] = td.measure_one(m, t, f_inj, tc, SET, all_mask, r_all=r_all)
                if any(v is None for v in results.values()):
                    continue
                i0 = min(np.searchsorted(t, tc), len(t) - 1)
                for name, r in results.items():
                    rows.append(dict(depth_inj=depth, shift=shift, n=n, tc=tc,
                                     since_orbit_start=float(since_start[i0]), method=name,
                                     recovered=r["depth_central"], recovered_shape=r["depth_shape"],
                                     bias=r["depth_central"] - depth, bias_shape=r["depth_shape"] - depth))
    return pd.DataFrame(rows)


def inj_summary(df):
    g = df.groupby(["background", "depth_inj", "method"])
    s = pd.DataFrame({
        "n_injections": g["bias"].count(),
        "mean_bias": g["bias"].mean(),
        "std_error_of_mean_bias": g["bias"].std() / np.sqrt(g["bias"].count()),
        "median_bias": g["bias"].median(),
        "scatter_std": g["bias"].std(),
        "rmse": g["bias"].apply(lambda b: float(np.sqrt(np.mean(b ** 2)))),
    }).reset_index()
    s["method"] = pd.Categorical(s["method"], order, ordered=True)
    return s.sort_values(["background", "depth_inj", "method"], ascending=[True, False, True])


# ------------------------------------------------------- real transits ----
real_tables, real_summaries, real_shape = [], [], []
for label, fl in flux_variants.items():
    df, unusable = real_depths(fl)
    df.insert(0, "background", label)
    real_tables.append(df)
    s = summarize_real(df, "depth_central").join(paired_vs_sap_none(df))
    s.insert(0, "background", label)
    real_summaries.append(s)
    sh = summarize_real(df, "depth_shape")
    sh.insert(0, "background", label)
    real_shape.append(sh)
real = pd.concat(real_tables)
real.to_csv(f"{outdir}/real_depths_per_transit.csv", index=False)
real_sum = pd.concat(real_summaries)
real_sum.to_csv(f"{outdir}/real_depth_summary.csv")
real_sum_shape = pd.concat(real_shape)
real_sum_shape.to_csv(f"{outdir}/real_depth_summary_shape_estimator.csv")

# --------------------------------------------- sensitivity on real data ----
sens_rows = []
for hw in (0.14, 0.17, 0.25, 0.30):
    for ch in (0.015, 0.02, 0.03):
        S2 = dict(SET, half_window=hw, central_half=ch)
        df, _ = real_depths(flux_variants["as observed"], settings=S2)
        for name, g in df.groupby("method"):
            sens_rows.append(dict(half_window=hw, central_half=ch, method=name,
                                  n=len(g), mean_depth=g.depth_central.mean()))
sens = pd.DataFrame(sens_rows)
sens.to_csv(f"{outdir}/real_sensitivity_window_and_central.csv", index=False)

# ------------------------------------------------ injection / recovery ----
INJ_DEPTHS = [0.0105, 0.003]
SHIFTS = [0.25, 0.75]
inj_tables = []
for label, fl in flux_variants.items():
    df = injections(fl, INJ_DEPTHS, SHIFTS)
    df.insert(0, "background", label)
    inj_tables.append(df)
inj = pd.concat(inj_tables)
inj.to_csv(f"{outdir}/injection_recovery_all.csv", index=False)
summary = inj_summary(inj)
summary.to_csv(f"{outdir}/injection_recovery_summary.csv", index=False)

inj["position"] = np.where(inj.since_orbit_start < 1.5, "first 1.5 d of orbit", "later in orbit")
g = inj.groupby(["background", "depth_inj", "position", "method"])["bias"]
bypos = pd.DataFrame({"n": g.count(), "mean_bias": g.mean(),
                      "rmse": g.apply(lambda b: float(np.sqrt(np.mean(b ** 2))))}).reset_index()
bypos.to_csv(f"{outdir}/injection_recovery_by_orbit_position.csv", index=False)

settings_out.update(
    n_real_transits_used=int(real.n.nunique()),
    real_transits_unusable=[[int(a), float(b)] for a, b in unusable],
    n_injections_per_depth=int(inj[(inj.depth_inj == INJ_DEPTHS[0]) & (inj.background == "as observed")]
                               .groupby(["shift", "n"]).ngroups),
    injection_depths=INJ_DEPTHS, injection_shifts=SHIFTS, methods=order,
    sap_over_pdcsap_median_ratio=float(np.median(d.sap_flux.values / d.pdcsap_flux.values)),
)
json.dump(settings_out, open(f"{outdir}/analysis_settings.json", "w"), indent=2)
np.save(f"{outdir}/phase_model.npy", phase_model)

pd.set_option("display.width", 220)
print(json.dumps({k: v for k, v in settings_out.items() if k != "methods"}, indent=1))
print("\nREAL DATA, central depth by method")
print(real_sum.round(6).to_string())
print("\nREAL DATA, shape estimator")
print(real_sum_shape.round(6).to_string())
print("\nINJECTION summary (depth 0.0105)")
print(summary[summary.depth_inj == 0.0105].round(6).to_string(index=False))
