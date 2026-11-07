# Data

## Source record

| Item | Value |
|---|---|
| Target | WASP-18 (planet WASP-18 b) |
| TIC ID | 100100827 |
| TESS sector | 2 |
| Pipeline and product | SPOC, light curve file (`_lc.fits`) |
| Exposure time | 120 seconds |
| Archive | MAST (Mikulski Archive for Space Telescopes) |
| Original file name | `tess2018234235059-s0002-0000000100100827-0121-s_lc.fits` |
| First search and download | 31 October 2026 |
| CSV export | 7 November 2026 |
| Search and download code | `lk.search_lightcurve("WASP-18", mission="TESS", author="SPOC", exptime=120)`, first result, then `.download(flux_column=...)` |
| Software in Colab | Lightkurve 2.6.0, NumPy 2.1.3, pandas 2.2.3, Matplotlib 3.10.0, Python 3.13 |
| Quality flags | Lightkurve default bitmask |

The search listed ten products for this target (sectors 2, 3, 29, 30, 69, 96 and 103 to 106). Only sector 2 is used.

## The CSV file

`wasp18_sector2.csv` has 18,299 rows and no missing values. It is not normalized.

| Column | Meaning |
|---|---|
| `time_btjd` | Time in days, BJD minus 2457000 |
| `sap_flux` | Simple aperture photometry flux (electrons per second). Not corrected for instrument trends |
| `sap_flux_err` | Reported uncertainty on `sap_flux` |
| `pdcsap_flux` | Pipeline corrected flux (electrons per second) |
| `pdcsap_flux_err` | Reported uncertainty on `pdcsap_flux` |

The cadence is 2 minutes. There is one gap of about 1.45 days between the two orbits of the sector (BTJD 1367.16 to 1368.61).

## Notes

* The FITS file is not stored here, because it can be downloaded again with the code above.
* The median of the PDCSAP flux is about 4 percent higher than the median of the SAP flux.
