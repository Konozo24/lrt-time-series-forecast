# 05_sarima.R - STANDALONE script.
# Topic: LRT Ampang line monthly ridership (trips/month), data.gov.my
# Daily Public Transport Ridership dataset (Prasarana Malaysia + Ministry
# of Transport, CC BY 4.0), 2019-01 to 2026-07 (91 months). SDG 11
# (Sustainable Cities and Communities), Target 11.2. Model family: ARIMA
# (seasonal).
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
# Role in the comparison: the seasonal member of the ARIMA family. EDA
# below confirms annual seasonality (STL seasonal component, clear
# ACF/PACF spike at lag 12), so seasonal terms are expected to pay for
# themselves here. The AICc gap against a non-seasonal ARIMA is the
# direct statistical evidence for that.
#
# NOT using auto.arima(): both the non-seasonal and seasonal orders are
# selected by an explicit grid search, so the search space and the
# selection criterion are visible rather than delegated to a stepwise
# black box.
#
# d and D are fixed FIRST by unit-root tests (KPSS for d, seasonal-
# strength test for D) and then held constant across the grid, because
# AICc is only comparable between models fitted to the same effective
# series - mixing different differencing orders into one ranking would
# be an invalid comparison. The selected orders are read dynamically off
# the grid (not hardcoded) since the search space here is 12 candidates,
# not a handful of binary flags.

# setup
pkgs <- c("forecast", "tseries", "dplyr", "lubridate", "ggplot2", "zoo", "Kendall", "arrow")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(arrow)      # read_parquet()
library(forecast)   # Arima(), accuracy(), checkresiduals(), ndiffs()/nsdiffs()
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

# Degrees of freedom consumed by the fitted ARMA parameters (p+q+P+Q,
# from fit$arma). Box.test()'s default fitdf = 0 treats residuals as if
# nothing had been estimated and makes the test too LENIENT - the
# p-value comes out larger than it should.
model_fitdf <- function(fit) {
  if (!is.null(fit$arma) && length(fit$arma) >= 4) sum(fit$arma[1:4]) else 0
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


# model: choose D (seasonal) then d (non-seasonal) by unit-root tests,
# then grid search p, q, P, Q with d and D fixed, ranked by AICc
D <- nsdiffs(train)
d <- ndiffs(if (D > 0) diff(train, lag = 12) else train, test = "kpss")
cat("Selected differencing: d =", d, ", D =", D, "(seasonal period 12)\n")

grid <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
grid_results <- data.frame()
for (i in seq_len(nrow(grid))) {
  fit <- tryCatch(
    Arima(train,
          order          = c(grid$p[i], d, grid$q[i]),
          seasonal       = list(order = c(grid$P[i], D, grid$Q[i]), period = 12),
          include.drift  = FALSE),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    grid_results <- rbind(grid_results, data.frame(
      p = grid$p[i], d = d, q = grid$q[i],
      P = grid$P[i], D = D, Q = grid$Q[i],
      AICc = fit$aicc, BIC = fit$bic))
  }
}
grid_results <- grid_results[order(grid_results$AICc), ]
cat("\n--- SARIMA candidate grid (top 8 by AICc) ---\n")
print(head(grid_results, 8), row.names = FALSE)

best  <- grid_results[1, ]
label <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[12]", best$p, best$d, best$q, best$P, best$D, best$Q)
cat("\nSelected:", label, " AICc =", round(best$AICc, 2), "\n")

fit_sarima <- Arima(train, order = c(best$p, best$d, best$q),
                    seasonal = list(order = c(best$P, best$D, best$Q), period = 12))
print(summary(fit_sarima))

fc_sarima <- forecast(fit_sarima, h = h)
acc_sarima <- accuracy(fc_sarima, test)
print(acc_sarima)


# residual diagnostics
resid_sarima <- residuals(fit_sarima)
fitdf_sarima <- model_fitdf(fit_sarima)

cat("\n== Ljung-Box on SARIMA residuals (want p > 0.05) ==\n")
lb12 <- Box.test(resid_sarima, lag = 12,      fitdf = fitdf_sarima, type = "Ljung-Box")
lb16 <- Box.test(resid_sarima, lag = LAG_MAX, fitdf = fitdf_sarima, type = "Ljung-Box")
print(lb12)
print(lb16)

cat("\n== Residual ACF lags outside the white-noise band ==\n")
n_out_12 <- acf_out_of_bounds(resid_sarima, lag.max = 12)
n_out_16 <- acf_out_of_bounds(resid_sarima, lag.max = LAG_MAX)
cat("n_lags_out_12:", n_out_12, " n_lags_out_16:", n_out_16, "\n")


# overfitting check: single holdout
mase_train <- acc_sarima["Training set", "MASE"]
mase_test  <- acc_sarima["Test set", "MASE"]
rmse_train <- acc_sarima["Training set", "RMSE"]
rmse_test  <- acc_sarima["Test set", "RMSE"]
mase_gap_holdout   <- abs(mase_test - mase_train) / mase_test
rmse_ratio_holdout <- rmse_test / rmse_train
cat("\nrmse_ratio_holdout:", round(rmse_ratio_holdout, 3),
    " mase_gap_pct_holdout:", round(mase_gap_holdout * 100, 1), "\n")


# overfitting check: rolling-origin CV (5 expanding-window folds, 12-month horizon each)
# Orders are re-selected per fold from a reduced grid (p,q = 0:2, P,Q =
# 0:1 - same as the main grid above), matching the shared pipeline's
# 09_rolling_cv.R.
n_full  <- length(y) # 91 months
origins <- seq(n_full - h - 16, n_full - h, by = 4)   # 5 expanding origins
cat("\nRolling-CV origins (training sizes):", paste(origins, collapse = ", "), "\n")

cv_sarima <- do.call(rbind, lapply(origins, function(ts_size) {
  tr <- head(y, ts_size)
  te <- window(y, start = time(y)[ts_size + 1], end = time(y)[ts_size + h])

  D_fold <- nsdiffs(tr)
  d_fold <- ndiffs(if (D_fold > 0) diff(tr, lag = 12) else tr, test = "kpss")
  g_fold <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
  r_fold <- data.frame()
  for (i in seq_len(nrow(g_fold))) {
    f <- tryCatch(Arima(tr, order = c(g_fold$p[i], d_fold, g_fold$q[i]),
                        seasonal = list(order = c(g_fold$P[i], D_fold, g_fold$Q[i]), period = 12)),
                  error = function(e) NULL)
    if (!is.null(f)) r_fold <- rbind(r_fold, data.frame(
      p = g_fold$p[i], q = g_fold$q[i], P = g_fold$P[i], Q = g_fold$Q[i], AICc = f$aicc))
  }
  r_fold <- r_fold[order(r_fold$AICc), ]
  m <- Arima(tr, order = c(r_fold$p[1], d_fold, r_fold$q[1]),
            seasonal = list(order = c(r_fold$P[1], D_fold, r_fold$Q[1]), period = 12))

  fcv <- forecast(m, h = h)
  a   <- accuracy(fcv, te)
  data.frame(train_size = ts_size,
             MAPE = a["Test set", "MAPE"],
             RMSE = a["Test set", "RMSE"],
             MASE = a["Test set", "MASE"])
}))
print(cv_sarima)

cv_summary <- cv_sarima %>%
  summarise(mean_MAPE = mean(MAPE), sd_MAPE = sd(MAPE),
            min_MAPE  = min(MAPE),  max_MAPE = max(MAPE),
            mean_RMSE = mean(RMSE), mean_MASE = mean(MASE),
            n_folds   = n())
print(cv_summary)

results <- data.frame(
  member        = "sarima",
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
write.csv(results, "output/member_sarima_results.csv", row.names = FALSE)

ggsave("output/plots/eda/stl_decomposition.png", p_stl,
       width = 9, height = 6, dpi = 150)
ggsave("output/plots/eda/mco_resolution.png", p_mco,
       width = 9, height = 5, dpi = 150)

# forecast vs actual, zoomed to the last 3 years
p_fc <- autoplot(fc_sarima) +
  autolayer(test, series = "Actual", color = "red") +
  ggtitle(paste(label, "- Forecast vs Actual")) +
  scale_y_millions() + xlim(2023, NA)
print(p_fc)
ggsave("output/plots/group_summary/fc_sarima.png", p_fc, width = 8, height = 5, dpi = 150)

png("output/plots/group_summary/resid_sarima.png", width = 900, height = 700, res = 130)
checkresiduals(fit_sarima)
dev.off()

cat("\nDone. Wrote output/member_sarima_results.csv and 3 plots to",
    "output/plots/group_summary/\n")
