# LRT Ampang Ridership Forecasting (BMMS2094)

Time series forecasting of LRT Ampang line monthly ridership.
SDG 11 (Sustainable Cities and Communities), Target 11.2.

Dataset: [Daily Public Transport Ridership](https://data.gov.my/data-catalogue/ridership_headline)
(data.gov.my, Prasarana Malaysia + Ministry of Transport, CC BY 4.0).

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
by AICc**, not by `auto.arima()` or the automatic `ets()` search — the
candidate sets and selection criteria are written out in each script and
printed to the console when it runs (`08_group_comparison.R`'s own
refits are also grid-searched the same way, independent of 04–07 - see
"Run order" below for why).

Differencing orders (`d`, `D`) are fixed first by unit-root tests
(`ndiffs`/`nsdiffs`) and then held constant across each grid. AICc is only
comparable between models fitted to the same effective series, so mixing
different differencing orders into one ranking would be invalid.

## Full period retained — MCO handling

All periods are included (2019-01 to 2026-07, 91 months); the series is
**not** truncated to remove the COVID-19 MCO disruption. The 2020–2021
collapse is resolved by **known-intervention imputation**: linear
interpolation of the trend across the disrupted window (March 2020 –
December 2021, defined from Malaysia's actual MCO timeline), plus a
seasonal component estimated by STL from the undisrupted months.

`02_mco_resolution.R` also documents an approach that was tried and
rejected first — automatic outlier detection by STL residual
thresholding — which produced false positives on genuinely normal months
adjacent to the shock, because a single sharp discontinuity distorts the
local trend estimate at its edges. Automatic detection suits scattered
outliers, not a sustained multi-month regime shift.

## Evaluation

- **Holdout**: last 12 months (one full seasonal cycle, so every calendar
  month appears exactly once in the test set).
- **Metrics**: MAPE, RMSE, MAE, MASE on the holdout.
- **Baseline**: SNAIVE — chosen over mean/naive/drift because the series
  has confirmed seasonality, so a benchmark that ignores seasonality
  would be an artificially weak comparison.
- **Residual diagnostics**: Ljung-Box at lag 24 (two seasonal cycles),
  plus a count of residual ACF lags exceeding ±1.96/√n.
- **Overfitting check**: train-vs-test accuracy gap must stay within 10%.
- **Robustness**: `09_rolling_cv.R` re-runs the comparison across 5
  rolling origins, since a single holdout is one draw.

## Run order

```
scripts/00_setup.R            # packages, shared settings and helpers
scripts/01_data_prep.R        # daily -> monthly, full period
scripts/02_mco_resolution.R   # resolve (not remove) the MCO window
scripts/03_eda.R              # decomposition, stationarity, ACF/PACF, lag plot
scripts/08_group_comparison.R # all four models + SNAIVE, one table
scripts/09_rolling_cv.R       # robustness check (slow — refits everything per fold)
```

That's the actual shared, dependent pipeline. Script 08 refits all four
models itself (ARIMA/SARIMA/ETS/BATS, each re-grid-searched inline) and
does not read anything from 04–07, so it only needs 02 to have run
first. Script 09 is independent of 04–08 the same way 08 is - it
re-grid-searches every model per fold.

`scripts/04_arima.R`, `05_sarima.R`, `06_ets.R`, and `07_bats.R` sit
**outside** this run order - despite the numbered filenames, none of
them `source("00_setup.R")` or depend on any other script's output.
Each is a fully self-contained, standalone version of one model (own
data read, own MCO resolution, own EDA, own diagnostics, own
holdout-and-rolling-CV overfitting checks), meant to be run and
submitted on its own as one member's individual work. Their selected
configs are the ones `08_group_comparison.R` and `09_rolling_cv.R`
hardcode/re-derive - see each script's own grid search for where that
choice comes from. Because they're independent of each other and of
08–10, running or skipping any one of them doesn't affect the rest of
the pipeline.

## Output structure

Two different sets of scripts write to the same `output/plots/` and
`output/tables/` folders, in two different styles, because 04–07 are
standalone (see "Run order" above) and 00–03/08–09 are the shared,
dependent pipeline.

**Shared pipeline (00–03, 08–09)** uses two small path helpers defined
in `00_setup.R` — `fig(category, "name.png")`, `tbl("name.csv")` —
instead of writing `"output/..."` paths directly:

```
output/
├── plots/
│   ├── eda/         EDA panels - decomposition, stationarity context
│   └── comparison/   cross-model plots (bar chart, CV stability,
│                      per-fold) - go in the GROUP report
└── tables/    model_comparison.csv, rolling_cv_*.csv, overfit_checks_cv.csv
               (all from 08/09)
```

Fitted model objects are not saved to disk in either style - each
script's `fit_*`/`fc_*` objects are used within that same run (plots,
diagnostics, the summary CSV) and nothing downstream reads a saved
model back in, so there is nothing to gain by persisting them.

**Standalone scripts (04–07)** write their own paths by hand (they
don't have `00_setup.R` loaded to provide the helpers), but land in the
same `output/plots/` and `output/tables/` folders as the shared pipeline:

```
output/
├── plots/
│   ├── eda/               stl_decomposition.png, mco_resolution.png -
│   │                       EACH of 04-07 writes to these SAME two
│   │                       filenames (see note below)
│   └── group_summary/     fc_<model>.png, resid_<model>.png - one pair
│                           per standalone script, distinct filenames
└── tables/
    ├── arima_result.csv    # one row: holdout + rolling-CV accuracy,
    ├── sarima_result.csv   #   overfitting checks, residual p-values
    ├── ets_result.csv      #   - written by 04/05/06/07 respectively
    └── bats_result.csv
```

Note on the shared filenames: all four of `04_arima.R`–`07_bats.R`
independently compute and save `plots/eda/stl_decomposition.png` and
`plots/eda/mco_resolution.png`. This is deliberate, not a bug - each
script reads the same source data and applies the same MCO
reconstruction logic, so the four scripts' versions should be
numerically identical; whichever one you run last is simply the one
left on disk. If you need more than one standalone script's copy to
coexist (e.g. comparing two members' output side by side), run them in
separate `output/` folders.

Key figures to pull into the report:
- `plots/eda/mco_resolution.png` — original vs. MCO-resolved series - the
  visual evidence for the whole reconstruction methodology. Produced by
  `02_mco_resolution.R` in the shared pipeline, or by whichever of
  `04_arima.R`–`07_bats.R` ran most recently if you're working from a
  standalone script instead.
- `plots/eda/eda_summary.png`, `stl_decomposition.png`, `seasonal_plot.png`,
  `lag_plot.png` — EDA (from `03_eda.R`; `stl_decomposition.png` is
  shared with 04–07 as noted above)
- `plots/group_summary/fc_<model>.png`, `resid_<model>.png` — one
  forecast + one residuals plot per model, for each member's individual
  report (from the standalone script for that model)
- `plots/comparison/model_comparison_plot.png` — all 4 models + SNAIVE vs.
  actual
- `plots/comparison/model_comparison_bar.png` — MAPE ranking bar chart
- `plots/comparison/rolling_cv_stability.png` — accuracy-vs-stability
  scatter (the ETS finding lives here: high mean, high variance)
- `plots/comparison/rolling_cv_per_fold.png` — MAPE per model across folds

## Importing into Posit Cloud from GitHub

1. Push to GitHub:
   ```
   cd "C:\Users\Ming\Desktop\Projects\lrt-time-series"
   git init
   git add .
   git commit -m "LRT Ampang forecasting: ARIMA/SARIMA/ETS/BATS"
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

Series: 91 months (2019-01 to 2026-07), train 79 / test 12.

Holdout accuracy, ranked by test MAPE (`output/tables/model_comparison.csv`):

| Model | MAPE | RMSE | MAE | MASE |
|---|---|---|---|---|
| ETS(A,N,A) | 2.898% | 246,123 | 178,010 | 0.198 |
| BATS | 3.302% | **226,856** | 196,828 | 0.219 |
| ARIMA(2,1,3) | 4.616% | 366,152 | 267,110 | 0.297 |
| SARIMA(0,1,2)(1,0,1)[12] | 4.711% | 317,058 | 278,865 | 0.310 |
| SNAIVE (baseline) | 5.218% | 396,504 | 322,437 | 0.358 |

Rolling-origin CV, 5 folds (`output/tables/rolling_cv_summary.csv`):

| Model | mean MAPE | sd | min | max | mean RMSE | mean MASE |
|---|---|---|---|---|---|---|
| **BATS** | **3.06%** | **0.81** | 2.33 | 4.35 | **215,748** | **0.202** |
| SARIMA | 4.08% | 0.50 | 3.38 | 4.71 | 274,150 | 0.266 |
| ETS | 5.23% | 2.39 | 2.90 | 8.01 | 364,209 | 0.346 |
| ARIMA | 5.90% | 2.45 | 4.06 | 10.07 | 411,745 | 0.390 |
| SNAIVE | 12.33% | 6.57 | 5.22 | 21.37 | 762,793 | 0.801 |

**Recommended model: BATS.** ETS wins the single holdout on MAPE, but is
third under cross-validation with a fold-to-fold spread three times BATS's
(sd 2.39 vs 0.81) - its holdout win does not survive being re-tested at
other origins. BATS is best on all three CV measures, wins RMSE on the
holdout itself, is never worse than 4.35% at any origin, and is the only
model to pass both overfitting rules against the CV mean (below).

All four models beat the SNAIVE benchmark (skill +0.10 to +0.45).

### Residual diagnostics

Ljung-Box at lag 12 and lag 16, with `fitdf = p+q+P+Q` so the p-values
match what `checkresiduals()` computes - `Box.test()`'s default of
`fitdf = 0` ignores the estimated parameters and is too lenient. Lag 16
rather than 24 because Hyndman's min(2m, T/5) rule caps at 15.8 for 79
training months.

| Model | fitdf | p (lag 12) | p (lag 16) | ACF out (12 / 16) | White noise? |
|---|---|---|---|---|---|
| SARIMA(0,1,2)(1,0,1)[12] | 4 | 0.412 | 0.388 | 1 / 1 | yes |
| BATS | 0 | 0.185 | 0.223 | 1 / 1 | yes |
| ETS(A,N,A) | 0 | 0.104 | 0.128 | 2 / 2 | yes |
| **ARIMA(2,1,3)** | 5 | **0.010** | **0.033** | 1 / 1 | **no** |

**ARIMA is the one model with significant residual autocorrelation**, and
it fails at both lags. A non-seasonal specification cannot capture all the
structure in monthly ridership - consistent with its higher AICc than
SARIMA (2230.80 vs 2215.56), and direct evidence that the seasonal terms
earn their place.

### Overfitting checks

Two rules: MASE gap <= 10%, and RMSE_test/RMSE_train <= 1.3. Both are
applied against the **CV mean** rather than the single holdout, since one
window is one draw (`output/tables/overfit_checks_cv.csv`).

| Model | RMSE ratio (CV) | <= 1.3 | MASE gap (CV) | <= 10% |
|---|---|---|---|---|
| **BATS** | **0.824** | yes | **0.5%** | **yes** |
| SARIMA | 0.863 | yes | 4.5% | yes |
| ARIMA | 1.204 | yes | 25.4% | no |
| ETS | 1.284 | yes | 30.3% | no |
| SNAIVE | 0.737 | yes | 24.8% | no |

Every model clears the 1.3 ratio, but ETS comes closest to breaching it
(1.284) and misses the MASE-gap rule by the widest margin of the four real
models. BATS is effectively flat between training and cross-validated
error (0.5% gap), i.e. it generalises without degradation.

On the single holdout the ratios are lower still - ETS 0.868, BATS 0.867,
SARIMA 0.998, ARIMA 1.071 - because train error there exceeds test error.
The CV column is the one to quote.

Training MAPE exceeds test MAPE for every model (e.g. ETS 4.88% vs
2.90%). That is not overfitting - the training window spans the MCO
collapse and the 2022-2024 recovery ramp, while the test window is a flat
mature period. It does mean the holdout figures are optimistic, which is
the reason the overfitting rules above are judged against the CV mean:
under cross-validation the direction reverses for ARIMA and ETS (ratios
1.204 and 1.284) and the holdout's flattering picture disappears.
