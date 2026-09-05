# 07_bats.R - STANDALONE script.
# Topic: LRT Ampang line monthly ridership (trips/month), data.gov.my
# Daily Public Transport Ridership dataset (Prasarana Malaysia + Ministry
# of Transport, CC BY 4.0), 2019-01 to 2026-07 (91 months). SDG 11
# (Sustainable Cities and Communities), Target 11.2. Model family: BATS
# (Box-Cox, ARMA errors, Trend, Seasonal - conventional/non-trigonometric
# seasonal representation, De Livera, Hyndman & Snyder, 2011).
#
# Self-contained: reads the raw dataset itself (via arrow::read_parquet
# straight off data.gov.my), resolves the COVID-19 MCO disruption itself,
# runs its own EDA and diagnostics, fits and validates its own model,
# writes its own outputs. No source(), no readRDS of a file another
# script produced, no dependency on any other script in this repo.
# Setup/data/MCO-resolution logic is duplicated from the shared pipeline
# on purpose so this one file can be submitted and run standalone.
#
# Why BATS rather than TBATS here: TBATS replaces the seasonal states
# with trigonometric (Fourier) terms, a representation designed for
# COMPLEX seasonality - multiple seasonal periods, non-integer periods,
# or very high frequency. This series has none of those: monthly
# ridership with a single integer period m = 12. With no complex
# seasonality to represent, BATS - which keeps the conventional
# seasonal-state representation - is the more direct specification.
# The trade-off, stated honestly: BATS estimates a full set of 12
# seasonal states where TBATS would use k <= 6 Fourier pairs, so BATS
# spends more parameters on a 79-observation training window. Whether
# that costs accuracy is exactly what the AIC comparison and holdout
# below are there to answer.
#
# Model pick: bats(use.trend = TRUE, use.box.cox = FALSE,
# use.damped.trend = TRUE, use.arma.errors = TRUE, seasonal.periods = 12),
# verified via an explicit 4-config grid search over {box-cox on/off} x
# {damped on/off}, ranked by AIC (not a single automatic bats() call).
#
#   BATS(no box-cox, damped)     AIC=2352.27, MAPE_test=3.30  <- picked
#   BATS(no box-cox, undamped)   AIC=2356.93, MAPE_test=6.37
#   BATS(box-cox, undamped)      AIC=2359.36, MAPE_test=5.97
#   BATS(box-cox, damped)        AIC=2363.02, MAPE_test=2.77
# AIC and holdout MAPE disagree here (box-cox+damped is more accurate
# out of sample), and that's reported below rather than hidden - but the
# AIC-best config is still what gets selected, because choosing by
# test-set accuracy would leak the holdout into model selection.


# setup
pkgs <- c("forecast", "tseries", "dplyr", "lubridate", "ggplot2", "zoo", "Kendall", "arrow")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(arrow)      # read_parquet()
library(forecast)   # bats(), accuracy(), checkresiduals()
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


# find the best AIC configuration {box-cox: on/off, damped: on/off}
configs <- expand.grid(bc = c(TRUE, FALSE), damped = c(TRUE, FALSE))

grid_results <- data.frame()
for (i in seq_len(nrow(configs))) {
  bc     <- configs$bc[i]
  damped <- configs$damped[i]
  fit <- tryCatch(
    bats(train, use.box.cox = bc, use.trend = TRUE, use.damped.trend = damped,
         use.arma.errors = TRUE, seasonal.periods = 12),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    fc <- forecast(fit, h = h)
    grid_results <- rbind(grid_results, data.frame(
      bc = bc, damped = damped,
      AIC       = fit$AIC,
      MAPE_test = accuracy(fc, test)["Test set", "MAPE"]))
  }
}
grid_results <- grid_results[order(grid_results$AIC), ]
cat("\n--- BATS candidate configurations (ranked by AIC) ---\n")
print(grid_results, row.names = FALSE)

# model
fit_bats <- bats(train, use.box.cox = FALSE, use.trend = TRUE, use.damped.trend = TRUE,
                  use.arma.errors = TRUE, seasonal.periods = 12)
print(fit_bats)
cat("Fitted Box-Cox lambda:", ifelse(is.null(fit_bats$lambda), NA, fit_bats$lambda), "\n")

fc_bats <- forecast(fit_bats, h = h)
acc_bats <- accuracy(fc_bats, test)
print(acc_bats)

# SNAIVE benchmark: any model that cannot beat SNAIVE is not earning the
# complexity it adds
snaive_fit <- snaive(train, h = h)
acc_snaive <- accuracy(snaive_fit, test)
cat("\n== SNAIVE benchmark (test-set accuracy) ==\n")
print(acc_snaive)
cat("BATS test MASE:", round(acc_bats["Test set", "MASE"], 3),
    "| SNAIVE test MASE:", round(acc_snaive["Test set", "MASE"], 3),
    "| BATS beats SNAIVE:", acc_bats["Test set", "MASE"] < acc_snaive["Test set", "MASE"], "\n")


# residual diagnostics
resid_bats <- residuals(fit_bats)

cat("\n== Ljung-Box on BATS residuals (want p > 0.05) ==\n")
lb12 <- Box.test(resid_bats, lag = 12,      type = "Ljung-Box")
lb16 <- Box.test(resid_bats, lag = LAG_MAX, type = "Ljung-Box")
print(lb12)
print(lb16)

cat("\n== Residual ACF lags outside the white-noise band ==\n")
n_out_12 <- acf_out_of_bounds(resid_bats, lag.max = 12)
n_out_16 <- acf_out_of_bounds(resid_bats, lag.max = LAG_MAX)
cat("n_lags_out_12:", n_out_12, " n_lags_out_16:", n_out_16, "\n")


# overfitting check: single holdout
mase_train <- acc_bats["Training set", "MASE"]
mase_test  <- acc_bats["Test set", "MASE"]
rmse_train <- acc_bats["Training set", "RMSE"]
rmse_test  <- acc_bats["Test set", "RMSE"]
mase_gap_holdout <- abs(mase_test - mase_train) / mase_test
rmse_ratio_holdout <- rmse_test / rmse_train
cat("\nrmse_ratio_holdout:", round(rmse_ratio_holdout, 3),
    " mase_gap_pct_holdout:", round(mase_gap_holdout * 100, 1), "\n")


# overfitting check: rolling-origin CV (5 expanding-window folds, 12-month horizon each)
n_full  <- length(y) # 91 months
origins <- seq(n_full - h - 16, n_full - h, by = 4)   # 5 expanding origins
cat("\nRolling-CV origins (training sizes):", paste(origins, collapse = ", "), "\n")

cv_bats <- do.call(rbind, lapply(origins, function(ts_size) {
  tr <- head(y, ts_size)
  te <- window(y, start = time(y)[ts_size + 1], end = time(y)[ts_size + h])
  m  <- bats(tr, use.box.cox = FALSE, use.trend = TRUE, use.damped.trend = TRUE,
             use.arma.errors = TRUE, seasonal.periods = 12)
  fcv <- forecast(m, h = h)
  a   <- accuracy(fcv, te)
  data.frame(train_size = ts_size,
             MAPE = a["Test set", "MAPE"],
             RMSE = a["Test set", "RMSE"],
             MASE = a["Test set", "MASE"])
}))
print(cv_bats)

cv_summary <- cv_bats %>%
  summarise(mean_MAPE = mean(MAPE), sd_MAPE = sd(MAPE),
            min_MAPE  = min(MAPE),  max_MAPE = max(MAPE),
            mean_RMSE = mean(RMSE), mean_MASE = mean(MASE),
            n_folds   = n())
print(cv_summary)

results <- data.frame(
  member        = "bats",
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
write.csv(results, "output/tables/bats_result.csv", row.names = FALSE)

ggsave("output/plots/eda/stl_decomposition.png", p_stl,
       width = 9, height = 6, dpi = 150)
ggsave("output/plots/eda/mco_resolution.png", p_mco,
       width = 9, height = 5, dpi = 150)

# forecast vs actual, zoomed to the last 3 years
p_fc <- autoplot(fc_bats) +
  autolayer(test, series = "Actual", color = "red") +
  ggtitle("BATS(no box-cox, damped) - Forecast vs Actual") +
  scale_y_millions() + xlim(2023, NA)
print(p_fc)
ggsave("output/plots/group_summary/fc_bats.png", p_fc, width = 8, height = 5, dpi = 150)

checkresiduals(fit_bats)   # shown on screen (Posit Cloud Plots pane)
dev.copy(png, filename = "output/plots/group_summary/resid_bats.png",
          width = 900, height = 700, res = 130)
dev.off()

cat("\nDone. Wrote output/tables/bats_result.csv, 2 plots to",
    "output/plots/eda/, and 2 plots to output/plots/group_summary/\n")
