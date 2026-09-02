# 05_sarima.R - Model family: ARIMA (seasonal).
#
# Role in the comparison: the seasonal member of the ARIMA family. EDA
# (03_eda.R) confirmed annual seasonality - STL seasonal strength ~0.48
# and a clear ACF/PACF spike at lag 12 - so seasonal terms are expected
# to pay for themselves here. The AICc gap against the non-seasonal
# ARIMA in 04_arima.R is the direct statistical evidence for that.
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
# be an invalid comparison.

source("scripts/00_setup.R")

y  <- load_series()
sp <- split_series(y)
train <- sp$train; test <- sp$test

# ---- Step 1: choose D (seasonal) then d (non-seasonal) -----------------
D <- nsdiffs(train)                                   # seasonal differences
d <- ndiffs(if (D > 0) diff(train, lag = 12) else train, test = "kpss")
cat("Selected differencing: d =", d, ", D =", D, "(seasonal period 12)\n")

# ---- Step 2: grid search p, q, P, Q with d and D fixed ----------------
grid <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
results <- data.frame()

for (i in seq_len(nrow(grid))) {
  fit <- tryCatch(
    Arima(train,
          order          = c(grid$p[i], d, grid$q[i]),
          seasonal       = list(order = c(grid$P[i], D, grid$Q[i]), period = 12),
          include.drift  = FALSE),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    results <- rbind(results, data.frame(
      p = grid$p[i], d = d, q = grid$q[i],
      P = grid$P[i], D = D, Q = grid$Q[i],
      AICc = fit$aicc, BIC = fit$bic))
  }
}

results <- results[order(results$AICc), ]
cat("\n--- SARIMA candidate grid (top 8 by AICc) ---\n")
print(head(results, 8), row.names = FALSE)
write.csv(results, tbl("sarima_grid.csv"), row.names = FALSE)

b <- results[1, ]
label <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[12]", b$p, b$d, b$q, b$P, b$D, b$Q)
cat("\nSelected:", label, " AICc =", round(b$AICc, 2), "\n")

# ---- Step 3: refit, forecast, evaluate --------------------------------
fit_sarima <- Arima(train, order = c(b$p, b$d, b$q),
                    seasonal = list(order = c(b$P, b$D, b$Q), period = 12))
print(summary(fit_sarima))

fc_sarima <- forecast(fit_sarima, h = H)
p_fc <- autoplot(fc_sarima) + autolayer(test, series = "Actual") +
  ggtitle(paste(label, "forecast"))
print(p_fc)
ggsave(fig("sarima_forecast.png"), p_fc, width = 8, height = 5, dpi = 150)

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_sarima, test))

# ---- Step 4: residual diagnostics -------------------------------------
cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_sarima)
dev.copy(png, filename = fig("sarima_residuals.png"), width = 900, height = 700, res = 130)
dev.off()
# Both lags: a spike at the seasonal lag stands out in a lag-12 test but
# gets diluted among mostly-zero terms in a longer one, so clear both.
# lb_test() supplies fitdf = p+q+P+Q, so these p-values agree with the
# "Model df" that checkresiduals() prints above; Box.test()'s default
# fitdf = 0 ignores the estimated parameters and is too lenient.
print(lb_test(fit_sarima, 12))
print(lb_test(fit_sarima, LAG_MAX))
cat("ACF-out-of-bounds: lag12 =", acf_out_of_bounds(residuals(fit_sarima), lag.max = 12),
    "/ 12 | lag", LAG_MAX, "=", acf_out_of_bounds(residuals(fit_sarima), lag.max = LAG_MAX),
    "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_sarima, test))
cat("\nOverfitting checks (", g$direction, ")\n")
cat("  MASE gap  :", round(g$mase_gap * 100, 1), "% | <= 10%  :", g$within_10pct, "\n")
cat("  RMSE ratio:", round(g$rmse_ratio, 3), "  | <= 1.3x :", g$within_1_3x, "\n")
cat("  (MAPE gap", round(g$mape_gap * 100, 1),
    "% reported only - distorted by the level difference between the\n",
    "  MCO-era training window and the test window, so not part of the rule)\n")

# ---- Step 5: save in the shared summary shape -------------------------
summ <- model_summary(label, fit_sarima, fc_sarima, test)
print(summ)
write.csv(summ, tbl("summary_sarima.csv"), row.names = FALSE)
saveRDS(fit_sarima, mdl("fit_sarima.rds"))
saveRDS(fc_sarima,  mdl("fc_sarima.rds"))
