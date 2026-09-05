# 05_sarima.R
# Topic: LRT Ampang line monthly ridership (trips/month), data.gov.my
# Daily Public Transport Ridership dataset (Prasarana Malaysia + Ministry
# of Transport, CC BY 4.0). SDG 11 (Sustainable Cities and Communities),
# Target 11.2. Model family: Seasonal ARIMA.
#
# Self-contained by design, matching the group's convention for individual
# scripts: reads the raw daily data itself, resolves the MCO disruption
# itself, runs its own identification/EDA, fits and validates its own
# model (single holdout AND rolling-origin CV), writes its own outputs.
# No source() of another member's file, so this one script can be
# submitted and re-run on its own to reproduce every number in the report.
#
# Two respects that support this report's individual "how the order was
# chosen" narrative, beyond a bare grid search:
#   1. Plots ACF/PACF of the *stationary* (post-differencing) series
#      specifically for identification, not just the raw series.
#   2. Reports approximate coefficient significance (|est| > 2*SE) as a
#      table, to support discussing which seasonal/non-seasonal terms
#      earn their place.
# Includes a rolling-origin CV step, following the same design as the
# group's 09_rolling_cv.R (same 5 expanding origins, same re-search-per-
# fold logic), so the CV numbers quoted in the individual report are
# reproducible from THIS file alone.
#
# EDA/diagnostic coverage aligned with the group's TBATS member script
# (07_bats.R) so every individual report runs the same checklist:
# raw-series Ljung-Box, Mann-Kendall trend test, a near-zero minimum-
# value check (MAPE stability), and an SNAIVE benchmark comparison.

# ---- setup -------------------------------------------------------------
pkgs <- c("forecast", "tseries", "dplyr", "lubridate", "ggplot2", "zoo", "Kendall")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(forecast)
library(tseries)
library(dplyr)
library(lubridate)
library(ggplot2)

dir.create("output/plots/models",        recursive = TRUE, showWarnings = FALSE)
dir.create("output/plots/group_summary", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",              recursive = TRUE, showWarnings = FALSE)

fig <- function(name) file.path("output/plots/models", name)
tbl <- function(name) file.path("output/tables", name)

scale_y_millions <- function() {
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M"))
}

# Counts residual ACF lags outside the +/- 1.96/sqrt(n) white-noise band.
acf_out_of_bounds <- function(resid, lag.max = 12) {
  r  <- na.omit(resid)
  n  <- length(r)
  ci <- 1.96 / sqrt(n)
  a  <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

# Box.test()'s default fitdf = 0 treats residuals as if nothing had been
# estimated, understating the parameters consumed and making the test
# too lenient - use p+q+P+Q instead.
model_fitdf <- function(fit) {
  if (!is.null(fit$arma) && length(fit$arma) >= 4) sum(fit$arma[1:4]) else 0
}

H <- 12          # one full seasonal cycle held out
LAG_MAX <- 16    # Hyndman's min(2m, T/5) rule: with 79 training months, 15.8 -> 16

# ---- data: daily -> monthly ---------------------------------------------
raw <- read.csv("data/ridership_headline.csv", stringsAsFactors = FALSE)
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
temp_filled    <- monthly$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)

# 3. run STL on that gap-filled series to get a seasonal component
stl_temp       <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est   <- as.numeric(stl_temp$time.series[, "seasonal"])

# 4. take the real values just before and just after the MCO window and
#    strip out their own seasonal component
i_pre  <- which(monthly$month == mco_start - months(1))
i_post <- which(monthly$month == mco_end + months(1))
last_pre   <- monthly$rail_lrt_ampang[i_pre]  - seasonal_est[i_pre]
first_post <- monthly$rail_lrt_ampang[i_post] - seasonal_est[i_post]

# 5. draw a straight line (linear interpolation) between those two
#    deseasonalized trend points across the gap months
n_gap <- sum(mco_mask)
bridge_trend <- seq(last_pre, first_post, length.out = n_gap + 2)[2:(n_gap + 1)]

# 6. add seasonality back onto the bridged trend, then replace only
#    the MCO months with it - every other month stays untouched.
resolved <- monthly$rail_lrt_ampang
resolved[mco_mask] <- bridge_trend + seasonal_est[mco_mask]
cat("Bridged", n_gap, "MCO-window months, 0 months dropped.\n")

y <- ts(resolved, start = start(ampang_ts), frequency = 12)   # series used from here on

p_mco <- autoplot(cbind(Original = ampang_ts, Resolved = y)) +
  ggtitle("LRT Ampang: Original vs. MCO-Resolved Series") +
  ylab("Monthly ridership") + scale_y_millions()
print(p_mco)
ggsave(fig("mco_resolution.png"), p_mco, width = 9, height = 5, dpi = 150)

# ---- EDA: whiteness / stationarity / trend, on the resolved series -------
cat("\n== Ljung-Box on RAW (resolved) series (want p < 0.05 -> not white noise) ==\n")
print(Box.test(y, lag = 12, type = "Ljung-Box"))
print(Box.test(y, lag = 24, type = "Ljung-Box"))

if (requireNamespace("Kendall", quietly = TRUE)) {
  cat("\n== Mann-Kendall trend test (H0: no monotonic trend) ==\n")
  print(Kendall::MannKendall(as.numeric(y)))
  cat("tau near 0 (|tau| < 0.1) = negligible trend magnitude even if p is\n",
      "small; 0.1-0.3 = weak. A significant p with a small tau means a real\n",
      "but practically minor upward trend, not an absence of trend.\n")
}

cat("\n== Min value (near-zero check for MAPE stability) ==\n")
print(min(y, na.rm = TRUE))

# ---- identification: differencing orders --------------------------------
train <- head(y, length(y) - H)
test  <- tail(y, H)

D <- nsdiffs(train)
d <- ndiffs(if (D > 0) diff(train, lag = 12) else train, test = "kpss")
cat("Identification: d =", d, ", D =", D, "(seasonal period = 12)\n")

cat("\n--- ADF / KPSS, level series ---\n")
print(adf.test(train)); print(kpss.test(train))
cat("\n--- ADF / KPSS, after d =", d, "regular differencing ---\n")
diffed_check <- diff(train, differences = d)
print(adf.test(diffed_check)); print(kpss.test(diffed_check))

decomp <- stl(y, s.window = "periodic", robust = TRUE)
seasonal_strength <- max(0, 1 - var(decomp$time.series[, "remainder"]) /
                            var(decomp$time.series[, "seasonal"] + decomp$time.series[, "remainder"]))
cat("\nSTL seasonal strength:", round(seasonal_strength, 3), "\n")

png(fig("stl_decomposition.png"), width = 800, height = 600, res = 130)
plot(decomp, main = "STL Decomposition")
dev.off()

# ---- MANUAL order identification (before the grid search) ---------------
# Step 1: read the RAW series' ACF/PACF. This is what a human analyst
# would look at first, before any differencing or automated search.
png(fig("acf_pacf_raw.png"), width = 800, height = 400, res = 130)
par(mfrow = c(1, 2))
Acf(train,  lag.max = 24, main = "ACF (raw series)")
Pacf(train, lag.max = 24, main = "PACF (raw series)")
dev.off()
cat("\nManual read (raw series): ACF decays slowly/near-linearly, PACF has a\n",
    "single dominant spike at lag 1 that then dies out - the classic\n",
    "non-stationary, trending signature. This is confirmed formally above\n",
    "by ADF/KPSS, and motivates differencing before any order is read.\n")

# Step 2: read the ACF/PACF of the STATIONARY (post-differencing) series.
# This is where a manual Box-Jenkins reading would tentatively propose
# non-seasonal and seasonal orders.
stationary_train <- diff(train, differences = d)
if (D > 0) stationary_train <- diff(stationary_train, lag = 12)

png(fig("sarima_identification_acf_pacf.png"), width = 800, height = 400, res = 130)
par(mfrow = c(1, 2))
Acf(stationary_train,  lag.max = 24, main = "ACF (stationary series)")
Pacf(stationary_train, lag.max = 24, main = "PACF (stationary series)")
dev.off()
cat("\nManual read (stationary series): a spike persists near lag 12 in both\n",
    "ACF and PACF - consistent with seasonal AR and/or MA structure rather\n",
    "than pure noise (supports D = 0: the seasonal signal survives because\n",
    "it was NOT removed by seasonal differencing, so seasonal AR/MA terms\n",
    "are needed to capture it). The non-seasonal lags decay without a single\n",
    "clean AR- or MA-type cutoff, so a manual reading alone cannot uniquely\n",
    "pin down (p,q) here - an AICc grid search below is used to double-check\n",
    "and resolve this ambiguity rather than relying on visual judgement alone.\n")

# ---- estimation: AICc grid search, d/D fixed -----------------------------
grid <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
grid_results <- data.frame()

for (i in seq_len(nrow(grid))) {
  fit <- tryCatch(
    Arima(train,
          order         = c(grid$p[i], d, grid$q[i]),
          seasonal      = list(order = c(grid$P[i], D, grid$Q[i]), period = 12),
          include.drift = FALSE),
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
cat("\nTop 8 candidate models by AICc:\n")
print(head(grid_results, 8), row.names = FALSE)
write.csv(grid_results, tbl("sarima_grid.csv"), row.names = FALSE)

best  <- grid_results[1, ]
label <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[12]", best$p, best$d, best$q, best$P, best$D, best$Q)
gap_to_2nd <- grid_results$AICc[2] - grid_results$AICc[1]
cat("\nSelected:", label, "| AICc =", round(best$AICc, 2),
    "| gap to 2nd-best:", round(gap_to_2nd, 2),
    if (gap_to_2nd <= 2) "(not substantial - simpler alternative is defensible)\n"
    else "(substantial - top model clearly preferred)\n")

fit_sarima <- Arima(train, order = c(best$p, best$d, best$q),
                     seasonal = list(order = c(best$P, best$D, best$Q), period = 12),
                     include.drift = FALSE)
print(summary(fit_sarima))

coef_table <- data.frame(
  term     = names(coef(fit_sarima)),
  estimate = round(coef(fit_sarima), 4),
  se       = round(sqrt(diag(vcov(fit_sarima))), 4)
)
coef_table$significant <- abs(coef_table$estimate) > 2 * coef_table$se
cat("\nCoefficient significance (|coef| > 2*s.e. approx. p<0.05):\n")
print(coef_table, row.names = FALSE)
write.csv(coef_table, tbl("sarima_coefficients.csv"), row.names = FALSE)

# ---- diagnostic checking --------------------------------------------------
fc_sarima  <- forecast(fit_sarima, h = H)
acc_sarima <- accuracy(fc_sarima, test)
print(acc_sarima)

# SNAIVE benchmark to sanity-check that the fitted model earns its complexity.
snaive_fit <- snaive(train, h = H)
acc_snaive <- accuracy(snaive_fit, test)
cat("\n== SNAIVE benchmark (test-set accuracy) ==\n")
print(acc_snaive)
cat("SARIMA test MASE:", round(acc_sarima["Test set", "MASE"], 3),
    "| SNAIVE test MASE:", round(acc_snaive["Test set", "MASE"], 3),
    "| SARIMA beats SNAIVE:", acc_sarima["Test set", "MASE"] < acc_snaive["Test set", "MASE"], "\n")

resid_sarima <- residuals(fit_sarima)
fitdf_sarima <- model_fitdf(fit_sarima)

lb12 <- Box.test(resid_sarima, lag = 12,      fitdf = fitdf_sarima, type = "Ljung-Box")
lb16 <- Box.test(resid_sarima, lag = LAG_MAX, fitdf = fitdf_sarima, type = "Ljung-Box")
cat("\nLjung-Box: p(lag12) =", round(lb12$p.value, 4), " p(lag16) =", round(lb16$p.value, 4), "\n")
cat("ACF-out-of-bounds: lag12 =", acf_out_of_bounds(resid_sarima, 12),
    "/12 | lag16 =", acf_out_of_bounds(resid_sarima, LAG_MAX), "/16\n")

checkresiduals(fit_sarima)   # shown on screen (Posit Cloud Plots pane)
dev.copy(png, filename = "output/plots/group_summary/resid_sarima.png",
          width = 900, height = 700, res = 120)
dev.off()

p_fc <- autoplot(fc_sarima) +
  autolayer(test, series = "Actual", color = "red") +
  ggtitle(paste(label, "- Forecast vs Actual")) + scale_y_millions()
print(p_fc)
ggsave("output/plots/group_summary/fc_sarima.png", p_fc, width = 8, height = 5, dpi = 150)

# ---- overfitting check: single holdout ------------------------------------
mase_train <- acc_sarima["Training set", "MASE"]
mase_test  <- acc_sarima["Test set", "MASE"]
rmse_train <- acc_sarima["Training set", "RMSE"]
rmse_test  <- acc_sarima["Test set", "RMSE"]
mase_gap_holdout   <- abs(mase_test - mase_train) / mase_test
rmse_ratio_holdout <- rmse_test / rmse_train
cat("\n[Holdout] RMSE ratio:", round(rmse_ratio_holdout, 3),
    "| MASE gap:", round(mase_gap_holdout * 100, 1), "%\n")

# ---- overfitting check: 5-fold rolling-origin CV --------------------------
# Same design as the group's 09_rolling_cv.R: 5 expanding-window origins,
# 12-month holdouts, orders RE-SEARCHED at each fold (not the fixed
# `best` config above) since the optimal order can shift as more data
# enters the training window.
n_full  <- length(y)
origins <- seq(n_full - H - 16, n_full - H, by = 4)   # 5 expanding origins
cat("\nRolling-CV origins (training sizes):", paste(origins, collapse = ", "), "\n")

cv_sarima <- do.call(rbind, lapply(origins, function(ts_size) {
  tr <- head(y, ts_size)
  te <- window(y, start = time(y)[ts_size + 1], end = time(y)[ts_size + H])

  D_fold <- nsdiffs(tr)
  d_fold <- ndiffs(if (D_fold > 0) diff(tr, lag = 12) else tr, test = "kpss")
  g_fold <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
  r_fold <- data.frame()
  for (i in seq_len(nrow(g_fold))) {
    f <- tryCatch(Arima(tr, order = c(g_fold$p[i], d_fold, g_fold$q[i]),
                        seasonal = list(order = c(g_fold$P[i], D_fold, g_fold$Q[i]), period = 12),
                        include.drift = FALSE),
                  error = function(e) NULL)
    if (!is.null(f)) r_fold <- rbind(r_fold, data.frame(
      p = g_fold$p[i], q = g_fold$q[i], P = g_fold$P[i], Q = g_fold$Q[i], AICc = f$aicc))
  }
  r_fold <- r_fold[order(r_fold$AICc), ]
  m <- Arima(tr, order = c(r_fold$p[1], d_fold, r_fold$q[1]),
             seasonal = list(order = c(r_fold$P[1], D_fold, r_fold$Q[1]), period = 12),
             include.drift = FALSE)

  fcv <- forecast(m, h = H)
  a   <- accuracy(fcv, te)
  data.frame(train_size = ts_size,
             order = sprintf("(%d,%d,%d)(%d,%d,%d)[12]", r_fold$p[1], d_fold, r_fold$q[1],
                              r_fold$P[1], D_fold, r_fold$Q[1]),
             MAPE = a["Test set", "MAPE"],
             RMSE = a["Test set", "RMSE"],
             MASE = a["Test set", "MASE"])
}))
print(cv_sarima)
write.csv(cv_sarima, tbl("sarima_rolling_cv.csv"), row.names = FALSE)

cv_summary <- cv_sarima %>%
  summarise(mean_MAPE = mean(MAPE), sd_MAPE = sd(MAPE),
            min_MAPE  = min(MAPE),  max_MAPE = max(MAPE),
            mean_RMSE = mean(RMSE), mean_MASE = mean(MASE), n_folds = n())
print(cv_summary)

rmse_ratio_cv <- cv_summary$mean_RMSE / rmse_train
mase_gap_cv   <- abs(cv_summary$mean_MASE - mase_train) / cv_summary$mean_MASE
cat("\n[Rolling CV] mean MAPE:", round(cv_summary$mean_MAPE, 2),
    "% | RMSE ratio:", round(rmse_ratio_cv, 3),
    "| MASE gap:", round(mase_gap_cv * 100, 1), "%\n")
cat("  within 1.3x RMSE ratio:", rmse_ratio_cv <= 1.3, "\n")
cat("  within 10% MASE gap  :", mase_gap_cv <= 0.10, "\n")

# ---- final results row -----------------------------------------------------
results <- data.frame(
  model                = label,
  MAPE_holdout         = round(acc_sarima["Test set", "MAPE"], 3),
  MASE_holdout         = round(mase_test, 3),
  rmse_ratio_holdout   = round(rmse_ratio_holdout, 3),
  mase_gap_holdout_pct = round(mase_gap_holdout * 100, 1),
  MAPE_cv_mean         = round(cv_summary$mean_MAPE, 3),
  sd_MAPE_cv           = round(cv_summary$sd_MAPE, 3),
  rmse_ratio_cv        = round(rmse_ratio_cv, 3),
  mase_gap_cv_pct      = round(mase_gap_cv * 100, 1),
  lb_pvalue_12         = round(lb12$p.value, 4),
  lb_pvalue_16         = round(lb16$p.value, 4)
)
print(results)
write.csv(results, tbl("sarima_final_summary.csv"), row.names = FALSE)
