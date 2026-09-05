# 04_arima.R - STANDALONE script.
# Topic: LRT Ampang line monthly ridership (trips/month), data.gov.my
# Daily Public Transport Ridership dataset (Prasarana Malaysia + Ministry
# of Transport, CC BY 4.0), 2019-01 to 2026-07 (91 months). SDG 11
# (Sustainable Cities and Communities), Target 11.2. Model family:
# STLARIMA (STL-decomposed ARIMA with intervention regressors).
#
# Self-contained: reads the raw dataset itself (via arrow::read_parquet
# straight off data.gov.my), resolves the COVID-19 MCO disruption itself,
# runs its own EDA and diagnostics, fits and validates its own model,
# writes its own outputs. No source(), no readRDS of a file another
# script produced, no dependency on any other script in this repo.
# Setup/data/MCO-resolution logic is duplicated from the shared pipeline
# on purpose so this one file can be submitted and run standalone.
#
# Role in the comparison: the ARIMA-family member that removes seasonality
# structurally (via STL) rather than modelling it with seasonal ARMA terms
# (SARIMA, 05_sarima.R) or exponential smoothing (ETS/BATS). The two xreg
# pulses give it a mechanism the other family members don't have: an
# explicit correction for one-off events at 2022-06 and 2023-04, rather
# than folding them into the noise.

# setup
pkgs <- c("forecast", "tseries", "dplyr", "lubridate", "ggplot2", "zoo", "Kendall", "arrow")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(arrow)      # read_parquet()
library(forecast)   # Arima(), forecast(), accuracy(), checkresiduals(), ndiffs(), seasadj()
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
# p-value comes out larger than it should. The xreg pulse coefficients
# are NOT counted here, matching the usual convention that fitdf refers
# to the ARMA part only (and matching checkresiduals()' own treatment).
model_fitdf <- function(fit) {
  if (!is.null(fit$arma) && length(fit$arma) >= 4) sum(fit$arma[1:4]) else 0
}

# Fits an ARIMA(p, d, q) + xreg model to a seasonally adjusted series and
# hands back a forecast/accuracy table on the ORIGINAL scale, by adding
# the STL seasonal component back onto the fitted values, the forecast
# mean, and the prediction interval bounds. accuracy()'s MASE needs the
# raw (not seasonally adjusted) training series in fc$x to scale correctly
# via seasonal differencing.
stlarima_accuracy_orig <- function(fit, fc, seasonal_train, seasonal_fc, tr_orig, te_orig) {
  fc$x      <- tr_orig
  fc$mean   <- ts(as.numeric(fc$mean) + seasonal_fc,
                   start = start(fc$mean), frequency = frequency(fc$mean))
  fc$fitted <- ts(as.numeric(fitted(fit)) + seasonal_train,
                   start = start(tr_orig), frequency = frequency(tr_orig))
  if (!is.null(fc$lower)) {
    fc$lower <- fc$lower + seasonal_fc
    fc$upper <- fc$upper + seasonal_fc
  }
  list(fc = fc, acc = accuracy(fc, te_orig))
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

# mco resolution: resolve the COVID-19 disruption without removing any
# months from the series (known-intervention imputation)

# 1. find the disrupted window - Malaysia's actual MCO timeline
mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")
mco_mask  <- monthly$month >= mco_start & monthly$month <= mco_end
cat("MCO window:", sum(mco_mask), "months (", format(mco_start), "to", format(mco_end), ")\n")

# 2. blank out that window and linearly fill the gap
temp_filled <- monthly$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)

# 3. run STL on that gap-filled series to get a seasonal component
stl_temp     <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est <- as.numeric(stl_temp$time.series[, "seasonal"])

# 4. take the real values just before and just after the MCO window and
#    strip out their own seasonal component
i_pre  <- which(monthly$month == mco_start - months(1))
i_post <- which(monthly$month == mco_end + months(1))
last_pre   <- monthly$rail_lrt_ampang[i_pre]  - seasonal_est[i_pre]     # seasonally adjusted
first_post <- monthly$rail_lrt_ampang[i_post] - seasonal_est[i_post]

# 5. draw a straight line (linear interpolation) between those two
#    deseasonalized trend points across the gap months
n_gap <- sum(mco_mask)
bridge_trend <- seq(last_pre, first_post, length.out = n_gap + 2)[2:(n_gap + 1)]

# 6. add seasonality back onto the bridged trend, then replace only
#    the MCO months with it - every other month stays untouched.
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

# train / test split - last 12 months (one full seasonal cycle)
h     <- 12
train <- head(y, length(y) - h)
test  <- tail(y, h)

# LAG_MAX for residual diagnostics: Hyndman's min(2m, T/5) rule. With
# T = 79 training months, 2m = 24 but T/5 = 15.8 binds first, so 16.
LAG_MAX <- 16


# EDA
p_stl <- autoplot(stl(y, s.window = "periodic", robust = TRUE)) +
  ggtitle("STL decomposition")
print(p_stl)

cat("\n== ADF on TRAIN (want p < 0.05 for stationary) ==\n")
print(adf.test(train))

cat("\n== KPSS on TRAIN (want p > 0.05 for stationary) ==\n")
print(kpss.test(train))

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


# STL-decompose the TRAINING series: model the seasonally adjusted part,
# carry the seasonal part forward untouched (repeat the last full year)
# to reconstruct forecasts on the original scale.
stl_train      <- stl(train, s.window = "periodic", robust = TRUE)
seasadj_train  <- seasadj(stl_train)
seasonal_train <- as.numeric(stl_train$time.series[, "seasonal"])
seasonal_fc_test <- rep(tail(seasonal_train, 12), length.out = h)

cat("\n== ADF on seas-adj train, level (want p < 0.05 = stationary) ==\n")
print(adf.test(seasadj_train))
cat("\n== KPSS on seas-adj train, level (want p > 0.05 = stationary) ==\n")
print(kpss.test(seasadj_train, null = "Level"))

diff_seasadj_train <- diff(seasadj_train)
cat("\n== ADF on seas-adj train, differenced ==\n")
print(adf.test(diff_seasadj_train))
cat("\n== KPSS on seas-adj train, differenced ==\n")
print(kpss.test(diff_seasadj_train, null = "Level"))

# ACF/PACF of the differenced, seasonally adjusted training series -
# the diagnostic that motivates the p/q search range below.
acf_vals  <- acf(diff_seasadj_train,  lag.max = 24, plot = FALSE)
pacf_vals <- pacf(diff_seasadj_train, lag.max = 24, plot = FALSE)
ci_bound  <- 1.96 / sqrt(length(diff_seasadj_train))

acf_df  <- data.frame(lag = round(as.numeric(acf_vals$lag) * 12), value = as.numeric(acf_vals$acf),
                       type = "ACF")
acf_df  <- acf_df[acf_df$lag > 0, ]
pacf_df <- data.frame(lag = round(as.numeric(pacf_vals$lag) * 12), value = as.numeric(pacf_vals$acf),
                       type = "PACF")
acf_pacf_df <- rbind(acf_df, pacf_df)
acf_pacf_df$type <- factor(acf_pacf_df$type, levels = c("ACF", "PACF"))

p_acf_pacf <- ggplot(acf_pacf_df, aes(x = lag, y = value)) +
  geom_col(fill = "#3366FF", width = 0.6) +
  geom_hline(yintercept = c(-ci_bound, ci_bound), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(title = "ACF/PACF - differenced, seasonally adjusted training series",
       x = "Lag (months)", y = "Correlation") +
  theme_minimal()
print(p_acf_pacf)


# known-event intervention regressors (xreg) - one-off pulse dummies at
# 2022-06 and 2023-04, per Justin's original specification.
intervention_dates <- as.Date(c("2022-06-01", "2023-04-01"))
xreg_all <- sapply(intervention_dates, function(d) as.numeric(monthly$month == d))
colnames(xreg_all) <- paste0("pulse_", format(intervention_dates, "%Y%m"))

xreg_train <- head(xreg_all, length(y) - h)
xreg_test  <- tail(xreg_all, h)
cat("\nIntervention months in training window:",
    paste(colnames(xreg_train)[colSums(xreg_train) > 0], collapse = ", "), "\n")


# model: choose d by unit-root test (KPSS) on the seasonally adjusted
# training series, then grid search p and q with d and xreg fixed,
# ranked by AICc on the training data only (no test-set leakage).
d <- ndiffs(seasadj_train, test = "kpss")
cat("Selected differencing order d =", d, "(KPSS-based, seas-adj series)\n")

grid <- expand.grid(p = 0:3, q = 0:3)
grid_results <- data.frame()
for (i in seq_len(nrow(grid))) {
  p <- grid$p[i]; q <- grid$q[i]
  fit <- tryCatch(
    Arima(seasadj_train, order = c(p, d, q), xreg = xreg_train),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    grid_results <- rbind(grid_results, data.frame(p = p, d = d, q = q,
                                                     AICc = fit$aicc, BIC = fit$bic))
  }
}
grid_results <- grid_results[order(grid_results$AICc), ]
cat("\n--- STLARIMA candidate grid (top 8 by AICc) ---\n")
print(head(grid_results, 8), row.names = FALSE)

best <- grid_results[1, ]
label <- paste0("STLARIMA(", best$p, ",", best$d, ",", best$q, ")+Interventions")
cat("\nSelected:", label, " AICc =", round(best$AICc, 2), "\n")

fit_arima <- Arima(seasadj_train, order = c(best$p, best$d, best$q),
                    xreg = xreg_train)
print(summary(fit_arima))

fc_arima_seasadj <- forecast(fit_arima, h = h, xreg = xreg_test)
recon <- stlarima_accuracy_orig(fit_arima, fc_arima_seasadj, seasonal_train, seasonal_fc_test, train, test)
fc_arima  <- recon$fc     # forecast object, reconstructed onto the original scale
acc_arima <- recon$acc
print(acc_arima)

# SNAIVE benchmark: any model that cannot beat SNAIVE is not earning the
# complexity it adds
snaive_fit <- snaive(train, h = h)
acc_snaive <- accuracy(snaive_fit, test)
cat("\n== SNAIVE benchmark (test-set accuracy) ==\n")
print(acc_snaive)
cat("ARIMA test MASE:", round(acc_arima["Test set", "MASE"], 3),
    "| SNAIVE test MASE:", round(acc_snaive["Test set", "MASE"], 3),
    "| ARIMA beats SNAIVE:", acc_arima["Test set", "MASE"] < acc_snaive["Test set", "MASE"], "\n")


# residual diagnostics (on the seasadj-scale model residuals, where the
# model was actually fit)
resid_arima <- residuals(fit_arima)
fitdf_arima <- model_fitdf(fit_arima)

cat("\n== Ljung-Box on ARIMA residuals (want p > 0.05) ==\n")
lb12 <- Box.test(resid_arima, lag = 12,      fitdf = fitdf_arima, type = "Ljung-Box")
lb16 <- Box.test(resid_arima, lag = LAG_MAX, fitdf = fitdf_arima, type = "Ljung-Box")
print(lb12)
print(lb16)

cat("\n== Residual ACF lags outside the white-noise band ==\n")
n_out_12 <- acf_out_of_bounds(resid_arima, lag.max = 12)
n_out_16 <- acf_out_of_bounds(resid_arima, lag.max = LAG_MAX)
cat("n_lags_out_12:", n_out_12, " n_lags_out_16:", n_out_16, "\n")


# overfitting check: single holdout
mase_train <- acc_arima["Training set", "MASE"]
mase_test  <- acc_arima["Test set", "MASE"]
rmse_train <- acc_arima["Training set", "RMSE"]
rmse_test  <- acc_arima["Test set", "RMSE"]
mase_gap_holdout   <- abs(mase_test - mase_train) / mase_test
rmse_ratio_holdout <- rmse_test / rmse_train
cat("\nrmse_ratio_holdout:", round(rmse_ratio_holdout, 3),
    " mase_gap_pct_holdout:", round(mase_gap_holdout * 100, 1), "\n")


# overfitting check: rolling-origin CV (5 expanding-window folds, 12-month
# horizon each). Order is re-selected per fold from a reduced grid
# (p,q = 0:2), matching the shared pipeline's 09_rolling_cv.R - a full
# p,q = 0:3 search at 5 origins would be needlessly slow for a
# robustness check. Each fold redoes its own STL decomposition and reuses
# the matching slice of xreg_all for that fold's train/test window.
n_full  <- length(y) # 91 months
origins <- seq(n_full - h - 16, n_full - h, by = 4)   # 5 expanding origins
cat("\nRolling-CV origins (training sizes):", paste(origins, collapse = ", "), "\n")

cv_arima <- do.call(rbind, lapply(origins, function(ts_size) {
  tr <- head(y, ts_size)
  te <- window(y, start = time(y)[ts_size + 1], end = time(y)[ts_size + h])

  xreg_tr <- xreg_all[1:ts_size, , drop = FALSE]
  xreg_te <- xreg_all[(ts_size + 1):(ts_size + h), , drop = FALSE]

  stl_fold      <- stl(tr, s.window = "periodic", robust = TRUE)
  seasadj_fold  <- seasadj(stl_fold)
  seasonal_fold <- as.numeric(stl_fold$time.series[, "seasonal"])
  seasonal_fc_fold <- rep(tail(seasonal_fold, 12), length.out = h)

  d_fold <- ndiffs(seasadj_fold, test = "kpss")
  g_fold <- expand.grid(p = 0:2, q = 0:2)
  r_fold <- data.frame()
  for (i in seq_len(nrow(g_fold))) {
    f <- tryCatch(Arima(seasadj_fold, order = c(g_fold$p[i], d_fold, g_fold$q[i]),
                         xreg = xreg_tr), error = function(e) NULL)
    if (!is.null(f)) r_fold <- rbind(r_fold, data.frame(p = g_fold$p[i], q = g_fold$q[i], AICc = f$aicc))
  }
  r_fold <- r_fold[order(r_fold$AICc), ]
  m <- Arima(seasadj_fold, order = c(r_fold$p[1], d_fold, r_fold$q[1]),
             xreg = xreg_tr)

  fcv <- forecast(m, h = h, xreg = xreg_te)
  a   <- stlarima_accuracy_orig(m, fcv, seasonal_fold, seasonal_fc_fold, tr, te)$acc
  data.frame(train_size = ts_size,
             MAPE = a["Test set", "MAPE"],
             RMSE = a["Test set", "RMSE"],
             MASE = a["Test set", "MASE"])
}))
print(cv_arima)

cv_summary <- cv_arima %>%
  summarise(mean_MAPE = mean(MAPE), sd_MAPE = sd(MAPE),
            min_MAPE  = min(MAPE),  max_MAPE = max(MAPE),
            mean_RMSE = mean(RMSE), mean_MASE = mean(MASE),
            n_folds   = n())
print(cv_summary)

results <- data.frame(
  member        = "arima",
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
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "output/tables/arima_result.csv", row.names = FALSE)
write.csv(grid_results, "output/tables/arima_grid.csv", row.names = FALSE)
write.csv(cv_arima, "output/tables/arima_rolling_cv.csv", row.names = FALSE)

ggsave("output/plots/eda/stl_decomposition.png", p_stl,
       width = 9, height = 6, dpi = 150)
ggsave("output/plots/eda/mco_resolution.png", p_mco,
       width = 9, height = 5, dpi = 150)
ggsave("output/plots/eda/acf_pacf_stlarima.png", p_acf_pacf,
       width = 8, height = 7, dpi = 150)

# forecast vs actual, zoomed to the last 3 years
p_fc <- autoplot(fc_arima) +
  autolayer(test, series = "Actual", color = "red") +
  ggtitle(paste(label, "- Forecast vs Actual")) +
  scale_y_millions() + xlim(2023, NA)
print(p_fc)
ggsave("output/plots/group_summary/fc_arima.png", p_fc, width = 8, height = 5, dpi = 150)

checkresiduals(fit_arima)   # shown on screen (Posit Cloud Plots pane)
dev.copy(png, filename = "output/plots/group_summary/resid_arima.png",
          width = 900, height = 700, res = 130)
dev.off()

cat("\nDone. Wrote 3 tables to output/tables/, 3 plots to output/plots/eda/, and 2 plots to",
    "output/plots/group_summary/\n")
