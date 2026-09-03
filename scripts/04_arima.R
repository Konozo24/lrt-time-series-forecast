# 04_arima.R - Model family: ARIMA (non-seasonal).
#
# Role in the comparison: the non-seasonal member of the ARIMA family.
# EDA (03_eda.R) checks for genuine annual seasonality (STL seasonal
# strength, ACF/PACF spike at lag 12) - see that script's console output
# for the current value. If seasonality is confirmed, a model with no
# seasonal terms is NOT expected to win. It is fitted deliberately as the
# within-family reference point: the AICc gap between this model and the
# seasonal SARIMA in 05_sarima.R is the direct evidence that the seasonal
# terms are earning their extra parameters, rather than being assumed.
#
# NOT using auto.arima(): orders are selected by an explicit grid search
# so the selection is visible and defensible, rather than delegated to a
# black-box stepwise search.
#
# Differencing order d is chosen FIRST, by unit-root test, and then held
# fixed across the grid. This matters: AICc is only comparable between
# models fitted to the same data, and changing d changes the effective
# series being modelled - so mixing different d values into one AICc
# ranking would be an invalid comparison.

source("scripts/00_setup.R")

y  <- load_series()
sp <- split_series(y)          # h = 12, one full seasonal cycle held out
train <- sp$train; test <- sp$test

# ---- Step 1: choose d by unit-root test (KPSS), not by guessing --------
d <- ndiffs(train, test = "kpss")
cat("Selected differencing order d =", d, "(KPSS-based)\n")

# ---- Step 2: grid search p and q with d fixed --------------------------
grid <- expand.grid(p = 0:3, q = 0:3)
results <- data.frame()

for (i in seq_len(nrow(grid))) {
  p <- grid$p[i]; q <- grid$q[i]
  fit <- tryCatch(
    Arima(train, order = c(p, d, q), include.drift = (d > 0)),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    results <- rbind(results, data.frame(p = p, d = d, q = q,
                                         AICc = fit$aicc, BIC = fit$bic))
  }
}

results <- results[order(results$AICc), ]
cat("\n--- ARIMA candidate grid (top 8 by AICc) ---\n")
print(head(results, 8), row.names = FALSE)
write.csv(results, tbl("arima_grid.csv"), row.names = FALSE)

best <- results[1, ]
cat("\nSelected: ARIMA(", best$p, ",", best$d, ",", best$q, ")  AICc =",
    round(best$AICc, 2), "\n")

# ---- Step 3: refit the selected order, forecast, evaluate --------------
fit_arima <- Arima(train, order = c(best$p, best$d, best$q),
                   include.drift = (best$d > 0))
print(summary(fit_arima))

fc_arima <- forecast(fit_arima, h = H)
p_fc <- autoplot(fc_arima) + autolayer(test, series = "Actual") +
  ggtitle(paste0("ARIMA(", best$p, ",", best$d, ",", best$q, ") forecast")) +
  scale_y_millions()
print(p_fc)
ggsave(fig("models", "arima_forecast.png"), p_fc, width = 8, height = 5, dpi = 150)

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_arima, test))

# ---- Step 4: residual diagnostics -------------------------------------
# checkresiduals() draws to whatever device is currently active, so this
# shows in the RStudio Plots pane as normal, then dev.copy() saves a copy
# of that same rendered plot to file - no need to call checkresiduals()
# twice (which would duplicate its console Ljung-Box text output).
cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_arima)
dev.copy(png, filename = fig("models", "arima_residuals.png"), width = 900, height = 700, res = 130)
dev.off()
# Both lags: a spike at the seasonal lag stands out in a lag-12 test but
# gets diluted among mostly-zero terms in a longer one, so clear both.
# lb_test() supplies fitdf = p+q+P+Q, so these p-values agree with the
# "Model df" that checkresiduals() prints above; Box.test()'s default
# fitdf = 0 ignores the estimated parameters and is too lenient.
print(lb_test(fit_arima, 12))
print(lb_test(fit_arima, LAG_MAX))
cat("ACF-out-of-bounds: lag12 =", acf_out_of_bounds(residuals(fit_arima), lag.max = 12),
    "/ 12 | lag", LAG_MAX, "=", acf_out_of_bounds(residuals(fit_arima), lag.max = LAG_MAX),
    "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_arima, test))
cat("\nOverfitting checks (", g$direction, ")\n")
cat("  MASE gap  :", round(g$mase_gap * 100, 1), "% | <= 10%  :", g$within_10pct, "\n")
cat("  RMSE ratio:", round(g$rmse_ratio, 3), "  | <= 1.3x :", g$within_1_3x, "\n")
cat("  (MAPE gap", round(g$mape_gap * 100, 1),
    "% reported only - not the pass/fail rule; see gap_check() in\n",
    "  00_setup.R for why)\n")

# ---- Step 5: save in the shared summary shape -------------------------
summ <- model_summary(paste0("ARIMA(", best$p, ",", best$d, ",", best$q, ")"),
                      fit_arima, fc_arima, test)
print(summ)
write.csv(summ, tbl("summary_arima.csv"), row.names = FALSE)
