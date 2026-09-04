# 06_ets.R - STANDALONE script.
# Topic: LRT Ampang line monthly ridership (trips/month), data.gov.my
# Daily Public Transport Ridership dataset (Prasarana Malaysia + Ministry
# of Transport, CC BY 4.0), 2019-01 to 2026-07 (91 months). SDG 11
# (Sustainable Cities and Communities), Target 11.2. Model family:
# Exponential Smoothing (ETS).
#
# Self-contained: reads the raw dataset itself (via arrow::read_parquet
# straight off data.gov.my), resolves the COVID-19 MCO disruption itself,
# runs its own EDA and diagnostics, fits and validates its own model,
# writes its own outputs. No source(), no readRDS of a file another
# script produced, no dependency on any other script in this repo.
# Setup/data/MCO-resolution logic is duplicated
# from the shared pipeline on purpose so this one file can be submitted
# and run standalone.
#
# Role in the comparison: the conventional exponential-smoothing member.
# Where ARIMA/SARIMA handle seasonality by differencing, ETS models it
# as explicit seasonal states updated by smoothing - a structurally
# different mechanism, which is why both families are worth comparing
# on the same series rather than assuming one representation suits the
# data.
#
# NOT using the automatic ets() model search: each candidate
# specification is fitted explicitly by name and ranked by AICc, so the
# comparison between additive and multiplicative seasonality (and
# between damped and undamped trend) is visible in the output rather
# than hidden inside an automatic selection step.
#
# ETS notation is ETS(Error, Trend, Seasonal):
#   Error    A = additive, M = multiplicative
#   Trend    N = none, A = additive, Ad = additive damped
#   Seasonal N = none, A = additive, M = multiplicative
# The additive-vs-multiplicative seasonal comparison is the substantive
# one here: multiplicative means the seasonal swing scales with the
# level of ridership rather than staying a fixed number of trips.

# setup
pkgs <- c("forecast", "tseries", "dplyr", "lubridate", "ggplot2", "zoo", "Kendall", "arrow")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(arrow)      # read_parquet()
library(forecast)   # ets(), accuracy(), checkresiduals()
library(tseries)    # adf.test(), kpss.test()
library(dplyr)
library(lubridate)
library(ggplot2)

# Counts residual ACF lags outside the +/- 1.96/sqrt(n) white-noise band.
acf_out_of_bounds <- function(resid, lag.max = 12) {
  r  <- na.omit(resid)
  n  <- length(r)
  ci <- 1.96 / sqrt(n)
  a  <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# format scientific notation (e.g. 2e+06) to millions, for plot axes
scale_y_millions <- function() {
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M"))
}


# download dataset and aggregate daily -> monthly.
raw <- as.data.frame(read_parquet("https://storage.data.gov.my/transportation/ridership_headline.parquet"))
raw$date <- as.Date(raw$date)

monthly <- raw %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(rail_lrt_ampang = sum(rail_lrt_ampang, na.rm = TRUE),
            n_days = n()) %>%
  ungroup() %>%
  filter(n_days >= 28)   # drop partial calendar months

cat("Monthly series:", nrow(monthly), "months,",
    format(min(monthly$month)), "to", format(max(monthly$month)), "\n")

start_year  <- year(min(monthly$month))
start_month <- month(min(monthly$month))
ampang_ts   <- ts(monthly$rail_lrt_ampang, start = c(start_year, start_month), frequency = 12)

# ---------------------------------------------------------------------
# MCO resolution: resolve the COVID-19 disruption WITHOUT removing any
# months from the series (known-intervention imputation, not automatic
# outlier detection - the latter was tried first and produced false
# positives on genuinely normal months adjacent to the shock, because a
# single sharp discontinuity distorts STL's local trend estimate right
# at its edges).
#
# 1. Explicitly define the disrupted window from Malaysia's actual MCO
#    timeline (domain knowledge, not statistical detection): March 2020
#    - December 2021 (MCO 1.0 through the National Recovery Plan /
#    endemic transition), 22 months.
# 2. Bridge the TREND across the window with linear interpolation
#    between the last stable pre-MCO value and first stable post-
#    recovery value, both seasonally adjusted first (raw anchors would
#    double-count those months' own seasonality - February is this
#    series' deepest seasonal trough, so anchoring on raw Feb 2020 would
#    drag the whole bridge down for no real reason).
# 3. Add back a SEASONAL component estimated by STL from the
#    undisrupted months only.
# 4. Replace only the 22 disrupted months with (bridged trend +
#    seasonal estimate). All months remain in the series - none
#    dropped - but the collapse-and-recovery no longer distorts
#    trend/seasonal estimation for the model fitted on this series.
mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")
mco_mask  <- monthly$month >= mco_start & monthly$month <= mco_end
cat("MCO window:", sum(mco_mask), "months (", format(mco_start), "to", format(mco_end), ")\n")

temp_filled <- monthly$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)
stl_temp     <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est <- as.numeric(stl_temp$time.series[, "seasonal"])

i_pre  <- which(monthly$month == mco_start - months(1))
i_post <- which(monthly$month == mco_end + months(1))
last_pre   <- monthly$rail_lrt_ampang[i_pre]  - seasonal_est[i_pre]     # seasonally adjusted
first_post <- monthly$rail_lrt_ampang[i_post] - seasonal_est[i_post]
n_gap <- sum(mco_mask)
bridge_trend <- seq(last_pre, first_post, length.out = n_gap + 2)[2:(n_gap + 1)]

resolved <- monthly$rail_lrt_ampang
resolved[mco_mask] <- bridge_trend + seasonal_est[mco_mask]

cat("Bridged", n_gap, "MCO-window months (trend interpolation + STL seasonal),",
    "0 months dropped.\n")

ampang_resolved <- ts(resolved, start = start(ampang_ts), frequency = 12)

p_mco <- autoplot(cbind(Original = ampang_ts, Resolved = ampang_resolved)) +
  ggtitle("LRT Ampang: Original vs. MCO-Resolved Series") +
  ylab("Monthly ridership") + scale_y_millions()
print(p_mco)

y <- ampang_resolved   # series used for everything from here on


# EDA
p_stl <- autoplot(stl(y, s.window = "periodic", robust = TRUE)) +
  ggtitle("STL decomposition")
print(p_stl)

cat("\n== ADF (want p < 0.05 for stationary) ==\n")
print(adf.test(y))

cat("\n== KPSS (want p > 0.05 for stationary) ==\n")
print(kpss.test(y))

cat("\n== Ljung-Box on RAW series (want p < 0.05 -> not white noise) ==\n")
print(Box.test(y, lag = 12, type = "Ljung-Box"))
print(Box.test(y, lag = 24, type = "Ljung-Box"))

if (requireNamespace("Kendall", quietly = TRUE)) {
  cat("\n== Mann-Kendall trend test (H0: no monotonic trend) ==\n")
  print(Kendall::MannKendall(as.numeric(y)))
  cat("tau near 0 (|tau| < 0.1) = negligible trend magnitude even if p is\n",
      "small; 0.1-0.3 = weak.\n")
}

cat("\n== Min value (near-zero check for MAPE stability) ==\n")
print(min(y, na.rm = TRUE))


# train / test split - last 12 months (one full seasonal cycle)
h     <- 12
train <- head(y, length(y) - h)
test  <- tail(y, h)

# LAG_MAX for residual diagnostics: Hyndman's min(2m, T/5) rule. With
# T = 79 training months, 2m = 24 but T/5 = 15.8 binds first, so 16.
LAG_MAX <- 16


# model: explicit candidate specifications, ranked by AICc. Combinations
# are restricted to numerically stable ones: forecast::ets() treats
# additive error with multiplicative seasonality (A,*,M) as unstable, so
# those are fitted with multiplicative error instead.
candidates <- list(
  list(label = "ETS(A,N,N)",  model = "ANN", damped = NULL),
  list(label = "ETS(A,A,N)",  model = "AAN", damped = FALSE),
  list(label = "ETS(A,Ad,N)", model = "AAN", damped = TRUE),
  list(label = "ETS(A,N,A)",  model = "ANA", damped = NULL),
  list(label = "ETS(A,A,A)",  model = "AAA", damped = FALSE),
  list(label = "ETS(A,Ad,A)", model = "AAA", damped = TRUE),
  list(label = "ETS(M,N,M)",  model = "MNM", damped = NULL),
  list(label = "ETS(M,A,M)",  model = "MAM", damped = FALSE),
  list(label = "ETS(M,Ad,M)", model = "MAM", damped = TRUE)
)

grid_results <- data.frame()
fits <- list()
for (cand in candidates) {
  fit <- tryCatch(
    ets(train, model = cand$model, damped = cand$damped),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    fits[[cand$label]] <- fit
    grid_results <- rbind(grid_results, data.frame(spec = cand$label,
                                                    AICc = fit$aicc, BIC = fit$bic))
  }
}
grid_results <- grid_results[order(grid_results$AICc), ]
cat("\n--- ETS candidate specifications (ranked by AICc) ---\n")
print(grid_results, row.names = FALSE)

best_label <- grid_results$spec[1]
cat("\nSelected:", best_label, " AICc =", round(grid_results$AICc[1], 2), "\n")

fit_ets <- fits[[best_label]]
print(summary(fit_ets))

fc_ets <- forecast(fit_ets, h = h)
acc_ets <- accuracy(fc_ets, test)
print(acc_ets)


# residual diagnostics
resid_ets <- residuals(fit_ets)

cat("\n== Ljung-Box on ETS residuals (want p > 0.05) ==\n")
# fitdf = 0 for ETS (the default), matching forecast::checkresiduals()'
# own treatment of exponential-smoothing fits.
lb12 <- Box.test(resid_ets, lag = 12,      type = "Ljung-Box")
lb16 <- Box.test(resid_ets, lag = LAG_MAX, type = "Ljung-Box")
print(lb12)
print(lb16)

cat("\n== Residual ACF lags outside the white-noise band ==\n")
n_out_12 <- acf_out_of_bounds(resid_ets, lag.max = 12)
n_out_16 <- acf_out_of_bounds(resid_ets, lag.max = LAG_MAX)
cat("n_lags_out_12:", n_out_12, " n_lags_out_16:", n_out_16, "\n")


# overfitting check: single holdout
mase_train <- acc_ets["Training set", "MASE"]
mase_test  <- acc_ets["Test set", "MASE"]
rmse_train <- acc_ets["Training set", "RMSE"]
rmse_test  <- acc_ets["Test set", "RMSE"]
mase_gap_holdout   <- abs(mase_test - mase_train) / mase_test
rmse_ratio_holdout <- rmse_test / rmse_train
cat("\nrmse_ratio_holdout:", round(rmse_ratio_holdout, 3),
    " mase_gap_pct_holdout:", round(mase_gap_holdout * 100, 1), "\n")


# overfitting check: rolling-origin CV (5 expanding-window folds, 12-month horizon each)
# Spec is re-selected per fold from a reduced 6-candidate list (drops the
# undamped-trend variants), matching the shared pipeline's 09_rolling_cv.R.
n_full  <- length(y) # 91 months
origins <- seq(n_full - h - 16, n_full - h, by = 4)   # 5 expanding origins
cat("\nRolling-CV origins (training sizes):", paste(origins, collapse = ", "), "\n")

cv_ets <- do.call(rbind, lapply(origins, function(ts_size) {
  tr <- head(y, ts_size)
  te <- window(y, start = time(y)[ts_size + 1], end = time(y)[ts_size + h])

  specs_fold <- list(c("ANN", NA), c("AAN", "TRUE"), c("ANA", NA),
                     c("AAA", "TRUE"), c("MNM", NA), c("MAM", "TRUE"))
  r_fold <- data.frame(); f_fold <- list()
  for (s in specs_fold) {
    dmp <- if (is.na(s[2])) NULL else as.logical(s[2])
    f <- tryCatch(ets(tr, model = s[1], damped = dmp), error = function(e) NULL)
    if (!is.null(f)) {
      k <- paste0(s[1], ifelse(is.na(s[2]), "", s[2])); f_fold[[k]] <- f
      r_fold <- rbind(r_fold, data.frame(key = k, AICc = f$aicc))
    }
  }
  r_fold <- r_fold[order(r_fold$AICc), ]
  m <- f_fold[[r_fold$key[1]]]

  fcv <- forecast(m, h = h)
  a   <- accuracy(fcv, te)
  data.frame(train_size = ts_size,
             MAPE = a["Test set", "MAPE"],
             RMSE = a["Test set", "RMSE"],
             MASE = a["Test set", "MASE"])
}))
print(cv_ets)

cv_summary <- cv_ets %>%
  summarise(mean_MAPE = mean(MAPE), sd_MAPE = sd(MAPE),
            min_MAPE  = min(MAPE),  max_MAPE = max(MAPE),
            mean_RMSE = mean(RMSE), mean_MASE = mean(MASE),
            n_folds   = n())
print(cv_summary)

results <- data.frame(
  member        = "ets",
  MASE_train    = mase_train,
  RMSE_train    = rmse_train,
  MASE_cv       = cv_summary$mean_MASE,
  RMSE_cv       = cv_summary$mean_RMSE,
  sd_MAPE_cv    = cv_summary$sd_MAPE,
  n_folds       = cv_summary$n_folds,
  rmse_ratio_cv = cv_summary$mean_RMSE / rmse_train,
  gap_pct_cv    = abs(cv_summary$mean_MASE - mase_train) / cv_summary$mean_MASE,
  lb_pvalue_12  = lb12$p.value,
  lb_pvalue_16  = lb16$p.value
)
results$within_1_3x_cv  <- results$rmse_ratio_cv <= 1.3
results$within_10pct_cv <- results$gap_pct_cv <= 0.10
print(results)


# outputs
dir.create("output/plots/eda", recursive = TRUE, showWarnings = FALSE)
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/member_ets_results.csv", row.names = FALSE)

ggsave("output/plots/eda/stl_decomposition.png", p_stl,
       width = 9, height = 6, dpi = 150)
ggsave("output/plots/eda/mco_resolution.png", p_mco,
       width = 9, height = 5, dpi = 150)

# forecast vs actual, zoomed to the last 3 years
p_fc <- autoplot(fc_ets) +
  autolayer(test, series = "Actual", color = "red") +
  ggtitle(paste(best_label, "- Forecast vs Actual")) +
  scale_y_millions() + xlim(2023, NA)
print(p_fc)
ggsave("output/plots/group_summary/fc_ets.png", p_fc, width = 8, height = 5, dpi = 150)

png("output/plots/group_summary/resid_ets.png", width = 900, height = 700, res = 130)
checkresiduals(fit_ets)
dev.off()

cat("\nDone. Wrote output/member_ets_results.csv and 3 plots to",
    "output/plots/group_summary/\n")
