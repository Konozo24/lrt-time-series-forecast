# 00_setup.R - packages. Run once per Posit Cloud session.
# Free tier: 1GB RAM. Do NOT install prophet (Stan compile fails/times`` out).

pkgs <- c("forecast", "tseries", "dplyr", "ggplot2", "lubridate", "gridExtra")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(forecast)   # bats(), accuracy(), checkresiduals()
library(tseries)    # adf.test(), kpss.test()
library(dplyr)
library(ggplot2)
library(lubridate)
library(gridExtra)  # grid.arrange() - combine multiple ggplot panels into one figure

# Ensure output/ (split into figures/tables/models) and data/ exist - git
# does not track empty directories
dir.create("data", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("output/models",  showWarnings = FALSE, recursive = TRUE)

# Path helpers so every script writes to the right subfolder without
# spelling out "output/figures/..." everywhere:
#   fig("x.png")    -> "output/figures/x.png"   (plots)
#   tbl("x.csv")    -> "output/tables/x.csv"    (grids, summaries, comparisons)
#   mdl("x.rds")    -> "output/models/x.rds"    (fitted model / forecast objects)
fig <- function(name) file.path("output/figures", name)
tbl <- function(name) file.path("output/tables",  name)
mdl <- function(name) file.path("output/models",  name)

# ---------------------------------------------------------------------
# SHARED SETTINGS - every model script uses these, so all four models are
# fitted and scored on exactly the same data and the same split.
# ---------------------------------------------------------------------

# h = 12: hold out one full seasonal cycle, so every calendar month is
# represented exactly once in the test set. A shorter window would make
# test accuracy depend on which half of the year happened to be held out.
H <- 12

# Ljung-Box / ACF lag. Hyndman's rule for seasonal data is
# min(2m, T/5): 2m = 24 would be two full seasonal cycles, but with
# T = 79 training months the T/5 cap binds at 15.8, so 16 it is. Using 24
# here would mean testing 24 autocorrelations from 79 observations, which
# dilutes any real signal across mostly-noise terms and costs the test
# power. 16 is also exactly what forecast::checkresiduals() picks for this
# series, so the summary tables and the console output agree.
# (A longer series would take 24 - e.g. at T = 540, T/5 = 108 and 2m binds
# instead. The cap is what differs here, not the rule.)
LAG_MAX <- 16

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

# Degrees of freedom consumed by the fitted parameters. Box.test() defaults
# to fitdf = 0, which treats residuals as if nothing had been estimated and
# makes the test too LENIENT - the p-value comes out larger than it should.
# For an ARIMA-class fit the ARMA orders count: fitdf = p + q + P + Q, taken
# from $arma (differencing and drift/mean do not count). This reproduces the
# "Model df" that forecast::checkresiduals() reports, so the summary tables
# and the console output agree instead of contradicting each other.
# ETS and BATS get 0, matching checkresiduals()' own treatment of them.
model_fitdf <- function(fit) {
  if (!is.null(fit$arma) && length(fit$arma) >= 4) sum(fit$arma[1:4]) else 0
}

# Ljung-Box with the correct degrees of freedom for this fit.
lb_test <- function(fit, lag) {
  Box.test(residuals(fit), lag = lag, fitdf = model_fitdf(fit),
           type = "Ljung-Box")
}

# Train-vs-test overfitting checks. TWO rules are reported:
#
#   1. MASE gap <= 10%   - the group's accuracy-gap convention.
#   2. RMSE ratio <= 1.3 - RMSE_test / RMSE_train. Unlike the gap this is SIGNED by
#                          construction: only test-worse-than-train can
#                          breach it, which is the direction overfitting
#                          actually points.
#
# The MAPE gap is still returned but is NOT part of the pass/fail rule.
# On this series it is badly distorted: the training window includes the
# MCO trough (~2.8m) while the test window sits near 6.1m, so dividing by
# a larger denominator shrinks test MAPE and inflates the gap for reasons
# that have nothing to do with overfitting. MASE is scale-free and RMSE
# is compared as a ratio, so both are immune to that.
#
# mape_gap/mase_gap use abs(), so they flag deviation in EITHER
# direction; `direction` records which way it actually went, because on
# this series train error usually EXCEEDS test error (the training window
# spans the disruption, the test window does not) - the opposite of
# overfitting, and not something a bare "FAIL" should hide.
gap_check <- function(acc_matrix) {
  mape_tr <- acc_matrix["Training set", "MAPE"]
  mape_te <- acc_matrix["Test set", "MAPE"]
  mase_tr <- acc_matrix["Training set", "MASE"]
  mase_te <- acc_matrix["Test set", "MASE"]
  rmse_tr <- acc_matrix["Training set", "RMSE"]
  rmse_te <- acc_matrix["Test set", "RMSE"]
  mape_gap   <- abs(mape_te - mape_tr) / mape_te
  mase_gap   <- abs(mase_te - mase_tr) / mase_te
  rmse_ratio <- rmse_te / rmse_tr
  list(mape_gap   = mape_gap,
       mase_gap   = mase_gap,
       rmse_ratio = rmse_ratio,
       rmse_train = rmse_tr,
       mase_train = mase_tr,
       direction  = if (mase_te > mase_tr) "test worse" else "train worse",
       within_10pct = (mase_gap <= 0.10),
       within_1_3x  = (rmse_ratio <= 1.3))
}

# One-line summary row per model, so all four scripts write results in an
# identical shape that 08_group_comparison.R can stack directly.
#
# Ljung-Box and the residual ACF are reported at BOTH lag 12 (one full
# seasonal cycle) and lag LAG_MAX = 16 (the min(2m, T/5) rule). One lag
# alone can mislead: a spike at the seasonal lag stands out clearly in a
# lag-12 test but gets diluted among mostly-zero terms in a longer one, so
# a model should clear both before its residuals are called white noise.
model_summary <- function(name, fit, fc, test) {
  a  <- accuracy(fc, test)
  g  <- gap_check(a)
  r  <- residuals(fit)
  lb12  <- lb_test(fit, 12)
  lbmax <- lb_test(fit, LAG_MAX)
  data.frame(
    model        = name,
    MAPE_test    = round(a["Test set", "MAPE"], 3),
    RMSE_test    = round(a["Test set", "RMSE"], 0),
    MAE_test     = round(a["Test set", "MAE"], 0),
    MASE_test    = round(a["Test set", "MASE"], 3),
    RMSE_train   = round(g$rmse_train, 0),
    MASE_train   = round(g$mase_train, 3),
    fitdf        = model_fitdf(fit),
    lb_pvalue_12 = round(lb12$p.value, 4),
    lb_pvalue_16 = round(lbmax$p.value, 4),
    lb_pass      = (lb12$p.value > 0.05) & (lbmax$p.value > 0.05),
    n_lags_out_12 = acf_out_of_bounds(r, lag.max = 12),
    n_lags_out_16 = acf_out_of_bounds(r, lag.max = LAG_MAX),
    mase_gap_pct = round(g$mase_gap * 100, 1),
    within_10pct = g$within_10pct,
    rmse_ratio   = round(g$rmse_ratio, 3),
    within_1_3x  = g$within_1_3x,
    direction    = g$direction,
    mape_gap_pct = round(g$mape_gap * 100, 1),   # reported, not a rule
    stringsAsFactors = FALSE
  )
}
