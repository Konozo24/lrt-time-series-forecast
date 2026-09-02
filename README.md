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
saved to `output/tables/*_grid.csv`.

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
scripts/04_arima.R
scripts/05_sarima.R
scripts/06_ets.R
scripts/07_bats.R
scripts/08_group_comparison.R # all four + SNAIVE, one table
scripts/09_rolling_cv.R       # robustness check (slow — refits everything per fold)
scripts/10_sensitivity_post_mco.R  # refit on 2022+ only (no reconstructed months)
```

Scripts 04–07 are independent of each other and can be run in any order
after 03. Script 08 refits everything itself, so it can also be run
standalone after 02. Script 10 reads `output/tables/model_comparison.csv`
for its side-by-side ranking, so run 08 before it.

## Output structure

Everything written by scripts 03–09 lands in one of three subfolders
under `output/`, so the folder doesn't fill up with 20+ mixed files:

```
output/
├── figures/   forecast plots, residual diagnostics, EDA panels, CV plots (.png)
├── tables/    grid searches, per-model summaries, comparisons (.csv)
└── models/    fitted model + forecast objects (.rds)
```

Each script uses three small path helpers defined in `00_setup.R` —
`fig("name.png")`, `tbl("name.csv")`, `mdl("name.rds")` — instead of
writing `"output/..."` paths directly, so the destination folder is
never spelled out by hand in each script.

Key figures to pull into the report:
- `figures/eda_summary.png`, `stl_decomposition.png`, `lag_plot.png` — EDA
- `figures/<model>_forecast.png`, `<model>_residuals.png` — per-model
- `figures/model_comparison_plot.png` — all 4 models + SNAIVE vs. actual
- `figures/model_comparison_bar.png` — MAPE ranking bar chart
- `figures/rolling_cv_stability.png` — accuracy-vs-stability scatter (the
  ETS finding lives here: high mean, high variance)
- `figures/rolling_cv_per_fold.png` — MAPE per model across folds

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

| Model | mean MAPE | sd | min | max |
|---|---|---|---|---|
| BATS | 3.06% | 0.81 | 2.33 | 4.35 |
| SARIMA | 4.08% | 0.50 | 3.38 | 4.71 |
| ETS | 5.23% | 2.39 | 2.90 | 8.01 |
| ARIMA | 5.90% | 2.45 | 4.06 | 10.07 |
| SNAIVE | 12.33% | 6.57 | 5.22 | 21.37 |

**Recommended model: BATS.** ETS wins the single holdout on MAPE, but is
third under cross-validation with a fold-to-fold spread three times BATS's
(sd 2.39 vs 0.81) - its holdout win does not survive being re-tested at
other origins. BATS is best on CV mean *and* wins RMSE on the holdout
itself, and it is the only model that is never worse than 4.35% at any
origin.

All four models beat the SNAIVE benchmark (skill +0.10 to +0.45).

Residual diagnostics use Ljung-Box at lag 12 and lag 16, with
`fitdf = p+q+P+Q` so the reported p-values match what
`checkresiduals()` computes. Lag 16 rather than 24 because Hyndman's
min(2m, T/5) rule caps at 15.8 for 79 training months:

| Model | fitdf | Ljung-Box p (lag 16) | White noise? |
|---|---|---|---|
| SARIMA(0,1,2)(1,0,1)[12] | 4 | 0.388 | yes |
| BATS | 0 | 0.223 | yes |
| ETS(A,N,A) | 0 | 0.128 | yes |
| **ARIMA(2,1,3)** | 5 | **0.033** | **no** |

**ARIMA is the one model with significant residual autocorrelation.** A
non-seasonal specification cannot capture all the structure in monthly
ridership - consistent with its higher AICc than SARIMA (2230.80 vs
2215.56), and direct evidence that the seasonal terms earn their place.

### Sensitivity: does the MCO reconstruction drive the result?

22 of the 79 training months are reconstructed rather than observed
(`02_mco_resolution.R`). `10_sensitivity_post_mco.R` refits the same
models on the 2022-01 onward subsample (55 months, 43 train / 12 test),
where no observation is reconstructed:

| Model | MAPE (post-MCO) | RMSE | MASE | Ljung-Box p | rank full → post |
|---|---|---|---|---|---|
| ETS | 4.168% | 318,488 | 0.294 | **0.0034** | 1 → 1 |
| BATS | 4.359% | 306,929 | 0.319 | 0.4760 | 2 → 2 |
| SARIMA | 5.042% | 355,867 | 0.358 | 0.3967 | 4 → 3 |
| SNAIVE | 5.218% | 396,504 | 0.391 | — | 5 → 4 |
| ARIMA | 9.971% | 676,968 | 0.715 | 0.3441 | 3 → 5 |

The top two are unchanged in both order and identity, so the headline
result does not depend on the reconstructed months. Two caveats: ARIMA
collapses on 43 observations without seasonal terms, and **ETS fails the
Ljung-Box test here (p = 0.0034)** - on unimputed data its accuracy edge
comes with autocorrelation left in the residuals, while BATS passes
cleanly (p = 0.476, 0 of 24 lags outside bounds). Rank correlation between
the two periods is 0.70. Note that 43 training months is only 3.6 seasonal
cycles, so this is a check on the *ranking*, not a fair accuracy contest.

Training MAPE exceeds test MAPE for every model (e.g. ETS 4.88% vs
2.90%). That is not overfitting - the training window spans the MCO
collapse and the 2022-2024 recovery ramp, while the test window is a flat
mature period. It does mean the holdout figures are optimistic.
