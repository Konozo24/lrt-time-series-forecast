# 00_setup.R - packages. Run once per Posit Cloud session.
# Free tier: 1GB RAM. Do NOT install prophet (Stan compile fails/times out).

pkgs <- c("forecast", "tseries", "dplyr", "ggplot2", "lubridate", "gridExtra")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(forecast)   # bats(), accuracy(), checkresiduals()
library(tseries)    # adf.test(), kpss.test()
library(dplyr)
library(ggplot2)
library(lubridate)
library(gridExtra)  # grid.arrange() - combine multiple ggplot panels into one figure

# Ensure output/ and data/ exist - git does not track empty directories.
# output/plots is split by AUDIENCE, matching how the figures get used when
# writing the report:
#   plots/eda/        EDA panels - decomposition, stationarity context
#   plots/comparison/ cross-model plots (bar chart, CV stability, per-fold) -
#                      these are what go in the GROUP report
# Each member's individual forecast/residuals plots come from their own
# standalone script (04-07), which write to plots/group_summary/ by hand.
dir.create("data", showWarnings = FALSE)
dir.create("output/plots/eda",        showWarnings = FALSE, recursive = TRUE)
dir.create("output/plots/comparison", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

# Path helpers so every script writes to the right subfolder without
# spelling out "output/plots/..." everywhere:
#   fig("eda", "x.png")   -> "output/plots/eda/x.png"
#   tbl("x.csv")          -> "output/tables/x.csv"  (grids, summaries, comparisons)
fig <- function(category, name) file.path("output/plots", category, name)
tbl <- function(name) file.path("output/tables", name)

# format scientific notation (e.g. 2e+06) to millions
scale_y_millions <- function() {
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M"))
}

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

# ---------------------------------------------------------------------
# STLARIMA HELPERS
# ---------------------------------------------------------------------
# 04_arima.R's model (credited to Justin, github.com/JusunF/ARIMA,
# combine.R): STL-decompose the series, fit ARIMA to the seasonally
# ADJUSTED series with two known-event pulse regressors (xreg), then add
# the seasonal component back to get forecasts on the original scale.
# 08_group_comparison.R and 09_rolling_cv.R share these helpers (rather
# than each hardcoding an order tuple the way they do for SARIMA/ETS/BATS)
# because a plain (p, d, q) triple isn't a complete spec here - it has to
# be paired with the STL decomposition and xreg it was selected against.
# See 04_arima.R for the full rationale, including why order selection is
# AICc-only (no test-set leakage).

INTERVENTION_DATES <- as.Date(c("2022-06-01", "2023-04-01"))

# Calendar date for every point of a monthly ts, derived purely from its
# start/frequency (the shared pipeline's ts objects don't carry the
# original data frame's month column around).
ts_month_dates <- function(y) {
  st <- start(y)
  seq(as.Date(sprintf("%d-%02d-01", st[1], st[2])), by = "month", length.out = length(y))
}

# One-off pulse dummy per intervention date, aligned to y.
intervention_xreg <- function(y) {
  dates <- ts_month_dates(y)
  m <- sapply(INTERVENTION_DATES, function(d) as.numeric(dates == d))
  colnames(m) <- paste0("pulse_", format(INTERVENTION_DATES, "%Y%m"))
  m
}

# STL-decompose tr, grid-search ARIMA(p,d,q)+xreg on the seasonally
# adjusted series (d chosen first by KPSS, p/q ranked by AICc), and
# return the fitted model plus what's needed to reconstruct forecasts on
# the original scale.
stlarima_fit <- function(tr, xreg_tr, p_range = 0:3, q_range = 0:3) {
  stl_tr      <- stl(tr, s.window = "periodic", robust = TRUE)
  seasadj_tr  <- seasadj(stl_tr)
  seasonal_tr <- as.numeric(stl_tr$time.series[, "seasonal"])
  d <- ndiffs(seasadj_tr, test = "kpss")

  grid <- expand.grid(p = p_range, q = q_range)
  res  <- data.frame()
  for (i in seq_len(nrow(grid))) {
    p <- grid$p[i]; q <- grid$q[i]
    f <- tryCatch(Arima(seasadj_tr, order = c(p, d, q), xreg = xreg_tr,
                        include.drift = (d > 0)), error = function(e) NULL)
    if (!is.null(f)) res <- rbind(res, data.frame(p = p, q = q, AICc = f$aicc))
  }
  res <- res[order(res$AICc), ]
  fit <- Arima(seasadj_tr, order = c(res$p[1], d, res$q[1]), xreg = xreg_tr,
               include.drift = (d > 0))
  list(fit = fit, seasonal_train = seasonal_tr, order = c(res$p[1], d, res$q[1]))
}

# Forecast an stlarima_fit() result h steps ahead and reconstruct it onto
# the ORIGINAL scale (add the seasonal component back to fitted/mean/PI),
# as a proper "forecast" object so accuracy()/autoplot() work on it
# directly, with MASE scaled correctly off the raw (not seasadj) series.
stlarima_forecast <- function(sf, tr, h, xreg_test) {
  fc <- forecast(sf$fit, h = h, xreg = xreg_test)
  seasonal_fc <- rep(tail(sf$seasonal_train, 12), length.out = h)
  fc$x      <- tr
  fc$mean   <- ts(as.numeric(fc$mean) + seasonal_fc,
                   start = start(fc$mean), frequency = frequency(fc$mean))
  fc$fitted <- ts(as.numeric(fitted(sf$fit)) + sf$seasonal_train,
                   start = start(tr), frequency = frequency(tr))
  if (!is.null(fc$lower)) {
    fc$lower <- fc$lower + seasonal_fc
    fc$upper <- fc$upper + seasonal_fc
  }
  fc
}
