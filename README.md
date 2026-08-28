# LRT Ampang Ridership Forecasting — TBATS (BMMS2094)

Individual contribution to BMMS2094 Statistics for Data Science group assignment.
SDG 11 (Sustainable Cities and Communities), Target 11.2.

Model family: ETS (TBATS variant).
Dataset: [Daily Public Transport Ridership](https://data.gov.my/data-catalogue/ridership_headline)
(data.gov.my, Prasarana Malaysia + Ministry of Transport, CC BY 4.0).

## Full period retained — MCO handling

Per tutor's instruction, **all periods are included** (2019-01 to 2026-06,
91 months) — the series is NOT truncated to exclude the COVID-19 MCO
disruption. Instead, the 2020-2021 collapse-and-recovery window is
resolved via a domain-knowledge-based intervention adjustment (linear
trend bridge + seasonal reconstruction over the known MCO window,
March 2020 – December 2021) rather than deleted. See
`scripts/02_mco_resolution.R` for the full method and reasoning, including
why a simpler automatic-outlier-detection approach was tried first and
rejected (it produced false positives on genuine pre/post-disruption
months).

## Run order

```
scripts/00_setup.R          # install/load packages
scripts/01_data_prep.R      # aggregate daily -> monthly, full period
scripts/02_mco_resolution.R # resolve (not remove) the MCO window
scripts/03_eda.R            # decomposition, stationarity, ACF/PACF
scripts/04_tbats_model.R    # fit TBATS, holdout accuracy, diagnostics
scripts/05_rolling_cv.R     # 7-fold rolling-origin CV
```

## Importing into Posit Cloud from GitHub

1. Push this folder to a new GitHub repo:
   ```
   cd "C:\Users\Ming\Desktop\Projects\lrt-time-series"
   git init
   git add .
   git commit -m "Initial TBATS project setup"
   git branch -M main
   git remote add origin https://github.com/<your-username>/lrt-time-series.git
   git push -u origin main
   ```
2. In Posit Cloud: **New Project → New Project from Git Repository**, paste
   the GitHub URL.
3. Open `lrt-time-series.Rproj` if it doesn't open automatically.
4. Run `scripts/00_setup.R` first each session (free tier — packages
   aren't persisted between long idle periods on some plans).
5. Run the remaining scripts in order (01 → 05).

## Reference numbers (from Python prototype — confirm R output lands close)

| Check | Expected result |
|---|---|
| Seasonal strength (STL) | ~0.48 |
| Best TBATS holdout MAPE (full, resolved series) | ~4.97% |
| Ljung-Box p-value (lag 6 / lag 12) | ~0.36 / ~0.38 (clean pass) |
| Rolling-CV mean MAPE | to be confirmed in R (Python check on truncated series: 5.17%) |

Numbers won't be identical between Python's `tbats` package and R's
`forecast::tbats()` (different optimizers) — small differences are
expected, large ones (order-of-magnitude) are not.
