"""Checks of the detrending tools on synthetic light curves with known answers."""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import transit_detrend as td

d = pd.read_csv(os.path.join(HERE, "..", "data", "wasp18_sector2.csv"))
t = d.time_btjd.values
P, T0, T14, T12 = 0.941459, 1354.45776, 0.0908, 0.0122
SET = dict(half_window=0.17, mask_half=0.75 * T14, central_half=0.02, t14=T14, t12=T12)
methods = td.make_methods()

def build(depth, ramp_amp, noise, seed=0, tau=0.7, eclipse=0.0):
    rng = np.random.default_rng(seed)
    f = np.ones_like(t)
    # ramp: starts at each orbit start (two orbits), exponential rise to 1, multiplicative
    starts = [t[0], t[np.argmax(np.diff(t) > 0.1) + 1]]
    since = np.where(t < starts[1], t - starts[0], t - starts[1])
    f *= 1.0 - ramp_amp * np.exp(-since / tau)
    x = td.time_from_event(t, P, T0)
    f *= 1.0 - depth * td.trapezoid_unit(x, T14, T12)
    xe = td.time_from_event(t, P, T0, shift=0.5)
    f *= 1.0 - eclipse * td.trapezoid_unit(xe, T14, T12)
    f += rng.normal(0, noise, t.size) if noise > 0 else 0
    return f

def event_masks():
    transit = td.event_mask(t, P, T0, SET["mask_half"])
    eclipse = td.event_mask(t, P, T0, SET["mask_half"], shift=0.5)
    return transit | eclipse

def run(f, tcs, mask):
    out = {}
    for m in methods:
        r_all = None
        if m["kind"] in ("none", "savgol"):
            r_all = td.detrend_global(m, t, f, mask)
        vals = []
        for tc in tcs:
            res = td.measure_one(m, t, f, tc, SET, mask, r_all=r_all)
            if res is not None:
                vals.append(res["depth_central"])
        out[m["name"]] = (np.mean(vals) if vals else np.nan, len(vals))
    return out

n_all = np.arange(-1, 31)
tcs = T0 + n_all * P
mask = event_masks()

print("TEST 1: no ramp, no noise, depth 0.0105. Every method should recover 0.0105")
res = run(build(0.0105, 0.0, 0.0), tcs, mask)
for k, (v, n) in res.items():
    print(f"  {k:22s} depth={v:.6f}  n={n}")
    assert n > 15, "too few usable transits"
    if "Savgol" not in k:
        assert abs(v - 0.0105) < 2e-5, k

print("\nTEST 2: ramp 0.004, no noise. 'none' must be biased, local poly degrees 2 and 3 close to truth")
res = run(build(0.0105, 0.004, 0.0), tcs, mask)
for k, (v, n) in res.items():
    print(f"  {k:22s} depth={v:.6f}  bias={v-0.0105:+.6f}")
assert abs(res["SAP none"][0] - 0.0105) > 5e-5

print("\nTEST 3: noise only (600 ppm). Scatter of the mean over 31 transits should be small")
vals = []
for seed in range(5):
    r = run(build(0.0105, 0.0, 6e-4, seed=seed), tcs, mask)
    vals.append(r["SAP local poly 1"][0])
print("  poly1 depth over 5 seeds:", np.round(vals, 5), "mean", round(np.mean(vals), 5))
assert abs(np.mean(vals) - 0.0105) < 1e-4

print("\nTEST 4: windows near the gap are rejected, not silently used")
gap_start = t[np.argmax(np.diff(t) > 0.1)]
mm = td.make_methods()[1]
res = td.measure_one(mm, t, build(0.0105, 0, 0), gap_start + 0.01, SET, mask)
print("  result for a window at the gap edge:", res)
assert res is None
res = td.measure_one(mm, t, build(0.0105, 0, 0), t[0] + 0.02, SET, mask)
print("  result for a window at the start of data:", res)
assert res is None

print("\nTEST 5: injection of a depth-0.003 transit at phase 0.25 is recovered when no ramp or noise")
f = build(0.0105, 0.0, 0.0)
tc = T0 + 10.25 * P
f_inj = td.inject_transit(t, f, tc, 0.003, T14, T12)
mask_inj = mask | td.event_mask(t, P, T0, SET["mask_half"], shift=0.25) & (np.abs(t - tc) < 0.5 * P)
m = np.zeros_like(mask); m |= mask; m |= np.abs(t - tc) < SET["mask_half"]
for meth in methods:
    r_all = td.detrend_global(meth, t, f_inj, m) if meth["kind"] in ("none", "savgol") else None
    res = td.measure_one(meth, t, f_inj, tc, SET, m, r_all=r_all)
    print(f"  {meth['name']:22s} recovered={res['depth_central']:.6f} shape={res['depth_shape']:.6f}")
    if "Savgol" not in meth["name"]:
        assert abs(res["depth_central"] - 0.003) < 5e-5

print("\nTEST 6: trapezoid shape properties")
x = np.linspace(-0.1, 0.1, 2001)
s = td.trapezoid_unit(x, T14, T12)
assert s.max() == 1.0 and s.min() == 0.0 and abs(s[1000] - 1.0) < 1e-12
assert abs((s > 0).sum() * (x[1] - x[0]) - T14) < 2e-3
print("  ok")

print("\nTEST 7: savgol on segments shorter than the window does not crash")
ts = np.arange(0, 0.2, td.CADENCE); fs = np.ones_like(ts)
r = td.savgol_flatten(ts, fs, np.zeros_like(ts, dtype=bool), 2.0)
assert np.allclose(r, 1.0)
print("  ok")
print("\nALL TESTS PASSED")
