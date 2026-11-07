"""Make the figures from the saved result tables. Run after run_analysis.py.

Usage: python src/make_figures.py <csv> <results folder> [figure folder]
"""
import sys, os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import transit_detrend as td

csv_path, outdir = sys.argv[1], sys.argv[2]
figdir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(outdir, "figures")
os.makedirs(figdir, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 10, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.labelsize": 10, "legend.fontsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
})
BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#7f7f7f"

S = json.load(open(f"{outdir}/analysis_settings.json"))
P, T0, T14 = S["period"], S["t0"], S["t14"]
d = pd.read_csv(csv_path)
t = d.time_btjd.values
sap = d.sap_flux.values / np.median(d.sap_flux.values)
pdc = d.pdcsap_flux.values / np.median(d.pdcsap_flux.values)
real = pd.read_csv(f"{outdir}/real_depths_per_transit.csv")
rsum = pd.read_csv(f"{outdir}/real_depth_summary.csv")
inj = pd.read_csv(f"{outdir}/injection_recovery_all.csv")
isum = pd.read_csv(f"{outdir}/injection_recovery_summary.csv")
order = S["methods"]
seg = td.split_segments(t)
starts = [t[a] for a, b in seg]

# ---- Figure 1: raw versus pipeline light curve and their ratio
fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
ax[0].plot(t, sap, ".", ms=1.2, color=BLUE, label="SAP (raw)")
ax[0].plot(t, pdc + 0.0, ".", ms=1.2, color=ORANGE, alpha=0.6, label="PDCSAP (pipeline corrected)")
ax[0].set_ylabel("Relative brightness")
ax[0].set_title("WASP-18, TESS sector 2")
ax[0].legend(loc="lower right", markerscale=6)
ratio = sap / pdc
ax[1].plot(t, (ratio / np.median(ratio) - 1) * 100, ".", ms=1.2, color=GREY)
ax[1].set_ylabel("SAP / PDCSAP (% change)")
ax[1].set_xlabel("Time (BTJD, days)")
for s0 in starts:
    for a in ax:
        a.axvline(s0, color=RED, lw=0.8, ls="--")
ax[1].text(starts[0] + 0.1, ax[1].get_ylim()[1] * 0.8, "orbit start", color=RED, fontsize=8)
fig.tight_layout(); fig.savefig(f"{figdir}/fig1_sap_vs_pdcsap.png"); plt.close(fig)

# ---- Figure 2: orbital brightness change (phase curve) check
real_mask = td.event_mask(t, P, T0, S["mask_half"]) | td.event_mask(t, P, T0, S["mask_half"], shift=0.5)
phase_model = np.load(f"{outdir}/phase_model.npy")
phi = ((t - T0) / P) % 1.0
bins = np.linspace(0, 1, 41)
idx = np.digitize(phi[~real_mask], bins) - 1
yb = np.array([np.mean(pdc[~real_mask][idx == k]) for k in range(40)])
eb = np.array([np.std(pdc[~real_mask][idx == k]) / np.sqrt((idx == k).sum()) for k in range(40)])
xb = 0.5 * (bins[1:] + bins[:-1])
order_ph = np.argsort(phi)
fig, ax = plt.subplots(figsize=(8, 4))
ax.errorbar(xb, (yb - 1) * 1e6, eb * 1e6, fmt="o", color=BLUE, ms=4, label="PDCSAP, binned, transit and eclipse excluded")
ax.plot(phi[order_ph], (phase_model[order_ph] - 1) * 1e6, color=RED, label="two harmonic fit")
ax.set_xlabel("Orbital phase (0 = transit)")
ax.set_ylabel("Brightness change (ppm)")
ax.set_title("Brightness change over one orbit")
ax.legend()
fig.tight_layout(); fig.savefig(f"{figdir}/fig2_orbital_brightness_change.png"); plt.close(fig)

# ---- Figure 3: depth of each real transit, three methods
fig, ax = plt.subplots(figsize=(9, 4.2))
obs = real[real.background == "as observed"]
for name, c, mk in [("SAP none", RED, "o"), ("SAP local poly 2", BLUE, "s"), ("PDCSAP none", ORANGE, "^")]:
    g = obs[obs.method == name]
    ax.plot(g.tc, g.depth_central * 100, mk + "-", color=c, ms=4, lw=0.8, label=name)
for s0 in starts:
    ax.axvline(s0, color=GREY, lw=0.8, ls="--")
ax.set_xlabel("Mid transit time (BTJD, days)")
ax.set_ylabel("Measured depth (%)")
ax.set_title("Depth of each transit by method")
ax.legend()
fig.tight_layout(); fig.savefig(f"{figdir}/fig3_depth_per_transit.png"); plt.close(fig)

# ---- Figure 4: mean real depth per method, with and without orbital brightness change
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True, sharex=True)
for ax, bg in zip(axes, ["as observed", "orbital brightness change removed"]):
    g = rsum[rsum.background == bg].set_index("method").loc[order]
    y = np.arange(len(order))[::-1]
    colors = [RED if m == "SAP none" else (ORANGE if m.startswith("PDCSAP") else BLUE) for m in order]
    ax.errorbar(g.mean_depth * 100, y, xerr=g.std_error * 100, fmt="none", ecolor=GREY, capsize=2)
    ax.scatter(g.mean_depth * 100, y, c=colors, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.set_xlabel("Mean central depth (%)")
    ax.set_title(bg.capitalize())
fig.suptitle("Real WASP-18 b depth by method (25 transits). Red: raw, blue: our detrending of SAP, orange: PDCSAP. Error bars: standard error", fontsize=9.5)
fig.tight_layout(); fig.savefig(f"{figdir}/fig4_real_depth_by_method.png"); plt.close(fig)

# ---- Figure 5: injection recovery, error per method
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
for ax, bg in zip(axes, ["as observed", "orbital brightness change removed"]):
    g = isum[(isum.background == bg) & (isum.depth_inj == 0.0105)].set_index("method").loc[order]
    y = np.arange(len(order))[::-1]
    ax.barh(y + 0.18, g.rmse * 1e6, height=0.34, color=BLUE, label="RMSE (total error)")
    ax.barh(y - 0.18, g.mean_bias.abs() * 1e6, height=0.34, color=ORANGE, label="|mean bias|")
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.set_xlabel("Error in recovered depth (ppm)")
    ax.set_title(bg.capitalize())
axes[0].legend(loc="lower right")
fig.suptitle("Injection recovery, injected depth 1.05 percent, 51 injections", fontsize=11)
fig.tight_layout(); fig.savefig(f"{figdir}/fig5_injection_error_by_method.png"); plt.close(fig)

# ---- Figure 6: injection bias versus time since orbit start
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for ax, bg in zip(axes, ["as observed", "orbital brightness change removed"]):
    for name, c, mk in [("SAP none", RED, "o"), ("SAP local poly 2", BLUE, "s"), ("SAP Savgol 1.0 d", GREEN, "^")]:
        g = inj[(inj.background == bg) & (inj.depth_inj == 0.0105) & (inj.method == name)]
        ax.scatter(g.since_orbit_start, g.bias * 1e6, c=c, marker=mk, s=14, label=name, alpha=0.85)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Days since orbit start")
    ax.set_title(bg.capitalize())
axes[0].set_ylabel("Recovered minus injected depth (ppm)")
axes[0].legend()
fig.suptitle("Where the errors come from", fontsize=11)
fig.tight_layout(); fig.savefig(f"{figdir}/fig6_injection_bias_vs_orbit_time.png"); plt.close(fig)

# ---- Figure 7: stacked transit after detrending
fig, ax = plt.subplots(figsize=(8, 4.2))
W = S["half_window"]
for name, c in [("SAP none", RED), ("SAP local poly 2", BLUE), ("PDCSAP none", ORANGE)]:
    m = [mm for mm in td.make_methods() if mm["name"] == name][0]
    fl = sap if m["flux"] == "sap" else pdc
    r_all = td.detrend_global(m, t, fl, real_mask) if m["kind"] in ("none", "savgol") else None
    xs, rs = [], []
    for n in real[real.background == "as observed"].n.unique():
        tc = T0 + n * P
        if m["kind"] == "local":
            out = td.local_poly_detrend(t, fl, tc, W, S["mask_half"], S["central_half"], m["degree"], extra_mask=real_mask)
            if out is None: continue
            win, r = out
            xs.append(t[win] - tc); rs.append(r)
        else:
            win = np.abs(t - tc) < W
            xs.append(t[win] - tc); rs.append(r_all[win])
    xs, rs = np.concatenate(xs), np.concatenate(rs)
    b = np.linspace(-W, W, 41)
    k = np.digitize(xs, b) - 1
    ax.plot(0.5 * (b[1:] + b[:-1]), [np.mean(rs[k == j]) for j in range(40)], "o-", ms=3, lw=0.8, color=c, label=name)
ax.axvline(-S["central_half"], color=GREY, ls=":"); ax.axvline(S["central_half"], color=GREY, ls=":")
ax.set_xlabel("Time from mid transit (days)")
ax.set_ylabel("Relative brightness")
ax.set_title("Stacked transit, 25 transits (dotted lines: depth measurement window)")
ax.legend()
fig.tight_layout(); fig.savefig(f"{figdir}/fig7_stacked_transit.png"); plt.close(fig)
print("figures written:", sorted(os.listdir(figdir)))
