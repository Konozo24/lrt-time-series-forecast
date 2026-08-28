# 00_setup.R - packages. Run once per Posit Cloud session.
# Free tier: 1GB RAM. Do NOT install prophet (Stan compile fails/times out).

pkgs <- c("forecast", "tseries", "dplyr", "ggplot2", "lubridate", "gridExtra")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(forecast)   # tbats(), accuracy(), checkresiduals()
library(tseries)    # adf.test(), kpss.test()
library(dplyr)
library(ggplot2)
library(lubridate)
library(gridExtra)  # grid.arrange() - combine multiple ggplot panels into one figure

# Ensure output/ and data/ exist - git does not track empty directories,
# so a fresh clone from GitHub won't have output/ until something creates
# it. Scripts that save plots/rds files into output/ will otherwise fail.
dir.create("output", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# Topic: LRT Ampang line monthly ridership, SDG 11 (Sustainable Cities and
# Communities), Target 11.2. Full period 2019-01 to 2026-06 retained per
# tutor's instruction (do NOT truncate out the MCO period) - see
# 02_mco_resolution.R for how the 2020-2021 disruption is handled instead
# of removed.

# ---------------------------------------------------------------------
# SHARED SETTINGS - every model script uses these, so all four models are
# fitted and scored on exactly the same data and the same split.
# ---------------------------------------------------------------------

# h = 12: hold out one full seasonal cycle, so every calendar month is
# represented exactly once in the test set. A shorter window would make
# test accuracy depend on which half of the year happened to be held out.
H <- 12

# Ljung-Box / ACF checks run to lag 24 = two full seasonal cycles for
# monthly data (Hyndman's rule of thumb for seasonal series: min(2m, n/5)).
LAG_MAX <- 24

load_series <- function() readRDS("data/ampang_monthly_full_resolved.rds")

split_series <- function(y, h = H) {
  list(train = head(y, length(y) - h), test = tail(y, h))
}

# ---------------------------------------------------------------------
# SHARED DIAGNOSTIC HELPERS
# ---------------------------------------------------------------------

# Count residual ACF lags (1..lag.max) exceeding +-1.96/sqrt(n). Roughly
# 5% of lags (~1 in 24) will exceed by chance even for genuine white
# noise; a much higher count means real structure was left unmodelled.
acf_out_of_bounds <- function(resid, lag.max = LAG_MAX) {
  r <- na.omit(resid)
  n <- length(r)
  ci <- 1.96 / sqrt(n)
  a <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# Train-vs-test accuracy gap, expressed as a fraction of the test value.
# Gap <= 10% on both MAPE and MASE is the group's agreed stability rule.
gap_check <- function(acc_matrix) {
  mape_tr <- acc_matrix["Training set", "MAPE"]
  mape_te <- acc_matrix["Test set", "MAPE"]
  mase_tr <- acc_matrix["Training set", "MASE"]
  mase_te <- acc_matrix["Test set", "MASE"]
  mape_gap <- abs(mape_te - mape_tr) / mape_te
  mase_gap <- abs(mase_te - mase_tr) / mase_te
  list(mape_gap = mape_gap, mase_gap = mase_gap,
       within_10pct = (mape_gap <= 0.10) & (mase_gap <= 0.10))
}

# One-line summary row per model, so all four scripts write results in an
# identical shape that 08_group_comparison.R can stack directly.
model_summary <- function(name, fit, fc, test) {
  a  <- accuracy(fc, test)
  g  <- gap_check(a)
  lb <- Box.test(residuals(fit), lag = LAG_MAX, type = "Ljung-Box")
  data.frame(
    model        = name,
    MAPE_test    = round(a["Test set", "MAPE"], 3),
    RMSE_test    = round(a["Test set", "RMSE"], 0),
    MAE_test     = round(a["Test set", "MAE"], 0),
    MASE_test    = round(a["Test set", "MASE"], 3),
    lb_pvalue    = round(lb$p.value, 4),
    n_lags_out   = acf_out_of_bounds(residuals(fit)),
    mape_gap_pct = round(g$mape_gap * 100, 1),
    within_10pct = g$within_10pct,
    stringsAsFactors = FALSE
  )
}
