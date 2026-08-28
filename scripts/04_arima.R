# 04_arima.R - Model family: ARIMA (non-seasonal).
#
# Role in the comparison: the non-seasonal member of the ARIMA family.
# EDA (03_eda.R) confirmed genuine annual seasonality (STL seasonal
# strength ~0.48, clear ACF/PACF spike at lag 12), so a model with no
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
write.csv(results, "output/arima_grid.csv", row.names = FALSE)

best <- results[1, ]
cat("\nSelected: ARIMA(", best$p, ",", best$d, ",", best$q, ")  AICc =",
    round(best$AICc, 2), "\n")

# ---- Step 3: refit the selected order, forecast, evaluate --------------
fit_arima <- Arima(train, order = c(best$p, best$d, best$q),
                   include.drift = (best$d > 0))
print(summary(fit_arima))

fc_arima <- forecast(fit_arima, h = H)
print(autoplot(fc_arima) + autolayer(test, series = "Actual") +
        ggtitle(paste0("ARIMA(", best$p, ",", best$d, ",", best$q, ") forecast")))

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_arima, test))

# ---- Step 4: residual diagnostics -------------------------------------
cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_arima)
print(Box.test(residuals(fit_arima), lag = LAG_MAX, type = "Ljung-Box"))
cat("ACF-out-of-bounds:", acf_out_of_bounds(residuals(fit_arima)), "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_arima, test))
cat("Train/test MAPE gap:", round(g$mape_gap * 100, 1),
    "% | within 10%:", g$within_10pct, "\n")

# ---- Step 5: save in the shared summary shape -------------------------
summ <- model_summary(paste0("ARIMA(", best$p, ",", best$d, ",", best$q, ")"),
                      fit_arima, fc_arima, test)
print(summ)
write.csv(summ, "output/summary_arima.csv", row.names = FALSE)
saveRDS(fit_arima, "output/fit_arima.rds")
saveRDS(fc_arima,  "output/fc_arima.rds")
