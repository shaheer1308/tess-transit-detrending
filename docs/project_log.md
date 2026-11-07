# Project Log

## Real data analysis settings

| Setting | Value |
|---|---|
| Data | WASP-18, TESS sector 2, SPOC, 120 s, 18,299 cadences (see `data/README.md`) |
| Normalization | Each flux column divided by its median |
| First period guess | Box Least Squares, 0.94155 days |
| Period (trapezoid fit to PDCSAP) | 0.941459 days |
| Mid transit time T0 | 1354.45776 BTJD |
| Total transit duration | 0.0908 days (2.18 hours) |
| Ingress time | 0.0122 days |
| Baseline mask half width | 0.0681 days (0.75 of the duration) |
| Local window half width | 0.17 days |
| Depth measuring half width | 0.02 days |
| Usable transits | 25 of 29 |
| Injection depths | 10,500 ppm and 3,000 ppm |
| Injection phases | 0.25 and 0.75 |
| Injections per depth | 51 |
| Random steps | None. Reruns are identical |

## Real data result: mean central depth over 25 transits

| Method | As observed (ppm) | Orbital brightness change removed (ppm) |
|---|---|---|
| SAP none | 10,854 | 10,443 |
| SAP local poly 0 | 10,596 | 10,356 |
| SAP local poly 1 | 10,596 | 10,356 |
| SAP local poly 2 | 10,370 | 10,327 |
| SAP local poly 3 | 10,371 | 10,327 |
| SAP Savgol 0.5 d | 10,417 | 10,320 |
| SAP Savgol 1.0 d | 10,628 | 10,337 |
| SAP Savgol 2.0 d | 10,770 | 10,361 |
| PDCSAP none | 10,922 | 10,504 |
| PDCSAP local poly 1 | 10,743 | 10,502 |

Standard errors, shape estimator results and window sensitivity are in `results/`.

## Injection recovery

See Table 3 of `docs/paper_draft.md` and `results/injection_recovery_summary.csv`.

## Checks against published numbers

* Our two harmonic fit of the orbital brightness change: cosine term -165 ppm, half period term -190 ppm, sine term 21 ppm. Published (Shporer et al. 2019): -174.4 ppm, -190.5 ppm and 21.0 ppm.
* Central depth expected from published parameters with limb darkening: 10,524 ppm. The squared radius ratio is 9,440 ppm. These are my own calculations from published numbers, not published depths.

## Earlier simulation results (run in Google Colab by the author)

Settings: 500 points, depth 0.05, noise 0.002, transit center 0.3 days, drift amplitude 0.01.

* Noise only, 100 trials: mean measured depth 0.05001, spread 0.00032.
* Drift added to flux, polynomial detrending: the mean error levelled off near 0.0004 below the true depth for degrees 3 to 8.
* Drift applied by multiplication: the levelling off disappeared. Degrees 4 to 8 gave errors between 0.00003 below and 0.00005 above the true depth.
* Best polynomial degree depended on drift period and on where the transit sat in the drift (Table 5 of the paper draft).

## Corrections made along the way

* An early message in the project said the dips in the light curve were about 1.2 days apart and gave a search range of 1.0 to 1.5 days. That was wrong. The true period is about 0.9415 days. A search over 0.5 to 3.0 days found it.
* A first Box Least Squares run on the narrow range returned 1.412 days. That is 1.5 times the true period, an alias, not a real signal.
* In the first injection results, low degree fits looked biased low by about 170 ppm. The cause is the planet's orbital brightness change, which has opposite curvature at the injection phases and at the transit. We therefore ran every injection test twice: once on the data as observed, and once after dividing out a fitted orbital brightness change.

## Open questions

* Why are PDCSAP depths about 1.5 percent higher than SAP depths after the orbital change is removed? One idea is the crowding correction. The value `CROWDSAP` in the FITS header can test this.
* Do the results hold in other sectors? Sector 3 is known to have stronger systematics (Shporer et al. 2019).
* How do the results change with a full transit model fit instead of the central mean depth?
* Why does the Savitzky Golay filter with a 0.5 day window show a small positive injection bias after the orbital change is removed?
