# Rapid KL Bus Ridership Forecasting (BMMS2094)

Time series forecasting of Rapid KL bus monthly ridership.
SDG 11 (Sustainable Cities and Communities), Target 11.2.

Dataset: [Daily Public Transport Ridership, bus_rkl column](https://data.gov.my/data-catalogue/ridership_headline?visual=bus_rkl)
(data.gov.my, Prasarana Malaysia + Ministry of Transport, CC BY 4.0).

**Coverage: 2022-01 to 2026-07 (55 months).** The `bus_rkl` column has no
data before 2022-01, so the series starts there - this is the column's
first observation, not a truncation of a longer one. Unlike an earlier
version of this project (LRT rail ridership, 2019 onward), this series
has no COVID-19 MCO-era gap to reconstruct: `bus_rkl` simply reports
nothing for the disrupted period, so every observation used here is a
real, unmodified measurement. There is no `02_mco_resolution.R` in this
version of the pipeline for that reason.

**Trade-off, stated honestly:** 55 months gives 43 training / 12 test
months - about 3.6 seasonal cycles, well short of the 6+ cycles usually
wanted for a seasonal model. Rolling-origin CV folds on this series go as
low as ~2.6 cycles, where SARIMA/BATS's seasonal terms have the least
room to estimate reliably (`09_rolling_cv.R`'s header documents the
narrower fold span chosen for that reason). Model comparisons here should
be read as a robustness check on which method degrades most gracefully
under a short sample, not a from-scratch accuracy benchmark.

## Models compared

Four models across two families, plus a benchmark:

| Script | Model | Family |
|---|---|---|
| `04_arima.R` | ARIMA(p,d,q) | ARIMA |
| `05_sarima.R` | SARIMA(p,d,q)(P,D,Q)[12] | ARIMA |
| `06_ets.R` | ETS(Error,Trend,Seasonal) | Exponential smoothing |
| `07_bats.R` | BATS | Exponential smoothing |
| `08_group_comparison.R` | SNAIVE | Benchmark (no estimated parameters) |

All orders and specifications are chosen by **explicit grid search ranked
by AICc** (AIC for BATS, which has no AICc), not by `auto.arima()` or the
automatic `ets()` search — the candidate sets and selection criteria are
written out in each script and saved to `output/tables/*_grid.csv`.

Differencing orders (`d`, `D`) are fixed first by unit-root tests
(`ndiffs`/`nsdiffs`) and then held constant across each grid. AICc is only
comparable between models fitted to the same effective series, so mixing
different differencing orders into one ranking would be invalid.

BATS's config (Box-Cox on/off × damped trend on/off) is re-selected live
in `07_bats.R` rather than hardcoded - see that script's header for why a
value carried over from a different dataset cannot be trusted here.
`08_group_comparison.R` re-selects it live too, for the same reason.
`09_rolling_cv.R` still uses a hardcoded placeholder pending the first run
of `07` - see the warning comment at that line.

## Evaluation

- **Holdout**: last 12 months (one full seasonal cycle, so every calendar
  month appears exactly once in the test set).
- **Metrics**: MAPE, RMSE, MAE, MASE on the holdout.
- **Baseline**: SNAIVE — chosen over mean/naive/drift because the series
  has confirmed seasonality, so a benchmark that ignores seasonality
  would be an artificially weak comparison.
- **Residual diagnostics**: Ljung-Box at lag 12 and `LAG_MAX` (currently 8
  - Hyndman's min(2m, T/5) rule for T = 43 training months; see the
  `LAG_MAX` comment in `00_setup.R`, and verify against what
  `checkresiduals()` reports on the first run), with `fitdf = p+q+P+Q` so
  the p-values match `checkresiduals()` rather than `Box.test()`'s
  default (too lenient). Plus a count of residual ACF lags exceeding
  ±1.96/√n at both lags.
- **Overfitting check**: two rules, both judged against the **rolling-CV
  mean** (the single holdout is one draw) - MASE gap ≤ 10%, and
  RMSE_test/RMSE_train ≤ 1.3.
- **Robustness**: `09_rolling_cv.R` re-runs the comparison across 5
  rolling origins, narrowed to a 12-month origin span (see that script's
  header) so the smallest fold still trains on ~2.6 seasonal cycles
  rather than ~2.25.

## Run order

```
scripts/00_setup.R            # packages, shared settings and helpers
scripts/01_data_prep.R        # daily -> monthly, bus_rkl, 2022-01 onward
scripts/03_eda.R              # decomposition, stationarity, ACF/PACF, lag plot
scripts/04_arima.R
scripts/05_sarima.R
scripts/06_ets.R
scripts/07_bats.R              # live BATS config grid - run before 08/09
scripts/08_group_comparison.R # all four + SNAIVE, one table
scripts/09_rolling_cv.R       # robustness check (slow — refits everything per fold)
```

(No `02` or `10`: this version of the pipeline has no MCO reconstruction
step and therefore nothing to sensitivity-check against it. Script
numbering keeps the gap rather than renumbering everything, so cross-
references in comments and old commits still resolve.)

Scripts 04–07 are independent of each other and can be run in any order
after 03. Script 08 refits everything itself, so it can also be run
standalone after 01. Run `07` before `08`/`09` so you know whether their
hardcoded/placeholder BATS config still needs updating.

## Output structure

Everything written by scripts 03–09 lands in one of two places under
`output/`. Plots are split by AUDIENCE - which report each one belongs
in - not just by script:

```
output/
├── plots/
│   ├── eda/         EDA panels - decomposition, stationarity context
│   ├── models/       one forecast + one residuals plot per model - go
│   │                  in each member's INDIVIDUAL report
│   └── comparison/   cross-model plots (bar chart, CV stability,
│                      per-fold) - go in the GROUP report
└── tables/    grid searches, per-model summaries, comparisons (.csv)
```

Each script uses two small path helpers defined in `00_setup.R` —
`fig(category, "name.png")`, `tbl("name.csv")` — instead of writing
`"output/..."` paths directly, so the destination folder is never
spelled out by hand in each script. Fitted model objects are not saved
to disk - each script's `fit_*`/`fc_*` objects are used within that
same run (plots, diagnostics, the summary CSV) and nothing downstream
reads a saved model back in, so there is nothing to gain by persisting
them.

Key figures to pull into the report:
- `plots/eda/eda_summary.png`, `stl_decomposition.png`, `seasonal_plot.png`,
  `lag_plot.png` — EDA
- `plots/models/<model>_forecast.png`, `<model>_residuals.png` — per-model,
  for individual reports
- `plots/comparison/model_comparison_plot.png` — all 4 models + SNAIVE vs.
  actual
- `plots/comparison/model_comparison_bar.png` — MAPE ranking bar chart
- `plots/comparison/rolling_cv_stability.png` — accuracy-vs-stability
  scatter
- `plots/comparison/rolling_cv_per_fold.png` — MAPE per model across folds

## Importing into Posit Cloud from GitHub

1. Push to GitHub:
   ```
   cd "C:\Users\Ming\Desktop\Projects\lrt-time-series"
   git init
   git add .
   git commit -m "Rapid KL bus forecasting: ARIMA/SARIMA/ETS/BATS"
   git branch -M main
   git remote add origin https://github.com/<your-username>/lrt-time-series.git
   git push -u origin main
   ```
2. Posit Cloud → **New Project → New Project from Git Repository**, paste
   the URL.
3. Run `scripts/00_setup.R` first each session.

Free tier is 1GB RAM. Do **not** install `prophet` — its Stan dependency
compiles from source and times out.

## Results

**Pending a fresh run.** The dataset changed from LRT rail ridership
(91 months, 2019-2026, MCO-reconstructed) to Rapid KL bus ridership
(55 months, 2022-2026, no reconstruction) - every number below has to be
regenerated on the new series before it can be reported. Run `00`→`09` on
Posit Cloud, then fill in this section from:

- `output/tables/model_comparison.csv` — holdout MAPE/RMSE/MAE/MASE per
  model, plus overfitting-rule and residual-diagnostic columns
- `output/tables/rolling_cv_summary.csv` — mean/sd/min/max MAPE, mean
  RMSE, mean MASE per model across the 5 CV folds
- `output/tables/overfit_checks_cv.csv` — the two overfitting rules
  judged against the CV mean (the version to quote, per Evaluation above)
- `output/tables/bats_grid.csv` — the 4-config BATS AIC/MAPE comparison
  from `07_bats.R`; check whether the AIC-best config matches what `08`
  and `09` are currently using and update them if not

Things to check specifically once real numbers exist, given the
short-sample trade-off noted above:

1. **Does any model breach the RMSE ratio ≤ 1.3 rule on the CV mean?**
   The equivalent post-2022 check on the old LRT series showed ARIMA and
   SARIMA breaching it at this sample length (ratios 2.42 and 1.33) -
   worth checking whether the same failure mode shows up here.
2. **Does any model lose to SNAIVE?** Same prior check showed ARIMA doing
   so (9.97% MAPE vs SNAIVE's 5.22%) on 43 training months.
3. **Do the rolling-CV folds actually fit?** With folds as short as ~2.6
   seasonal cycles, watch the console for `tryCatch` failures silently
   dropping candidates from a fold's grid search - a much smaller
   candidate pool for SARIMA/BATS in the earliest folds would be visible
   in `rolling_cv_folds.csv` as an unusually good-looking (under-searched)
   AICc winner.
4. **Confirm `LAG_MAX = 8` is still correct** by checking what
   `checkresiduals()` prints as "Total lags used" in each model's
   console output, per the comment in `00_setup.R`.
