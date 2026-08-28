# 07_tbats.R - Model family: Exponential Smoothing (TBATS variant).
#
# TBATS = Trigonometric seasonality, Box-Cox transformation, ARMA errors,
# Trend, Seasonal components (De Livera, Hyndman & Snyder, 2011).
#
# Role in the comparison: an extended exponential-smoothing state-space
# model, tested against the conventional ETS in 06_ets.R. The substantive
# difference is how seasonality is represented - ETS estimates 12 separate
# monthly seasonal states, while TBATS uses a small number of Fourier
# (trigonometric) terms, which is more parameter-efficient on a series
# this short. TBATS additionally offers a Box-Cox transform (to stabilise
# a seasonal swing that scales with the level) and an ARMA correction on
# the residuals; whether either is actually needed here is tested rather
# than assumed.
#
# NOT relying on a single automatic fit: the Box-Cox and damped-trend
# components are switched on and off explicitly across four candidate
# configurations, so their contribution is visible in the output.

source("scripts/00_setup.R")

y  <- load_series()
sp <- split_series(y)
train <- sp$train; test <- sp$test

# ---- Step 1: explicit component configurations ------------------------
configs <- list(
  list(label = "TBATS(box-cox, damped)",    bc = TRUE,  damped = TRUE),
  list(label = "TBATS(box-cox, undamped)",  bc = TRUE,  damped = FALSE),
  list(label = "TBATS(no box-cox, damped)", bc = FALSE, damped = TRUE),
  list(label = "TBATS(no box-cox, undamped)", bc = FALSE, damped = FALSE)
)

results <- data.frame()
fits <- list()

for (cfg in configs) {
  fit <- tryCatch(
    tbats(train, use.box.cox = cfg$bc, use.trend = TRUE,
          use.damped.trend = cfg$damped),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    fits[[cfg$label]] <- fit
    fc  <- forecast(fit, h = H)
    results <- rbind(results, data.frame(
      config    = cfg$label,
      AIC       = fit$AIC,
      MAPE_test = accuracy(fc, test)["Test set", "MAPE"]))
  }
}

results <- results[order(results$AIC), ]
cat("\n--- TBATS candidate configurations ---\n")
print(results, row.names = FALSE)
write.csv(results, "output/tbats_grid.csv", row.names = FALSE)

# Selection note: AIC differences between these configurations are small,
# so the damped/undamped choice is decided on out-of-sample accuracy
# rather than in-sample fit alone - damping constrains how far the trend
# is extrapolated over a 12-month horizon, which matters more for
# forecasting than the marginal AIC difference suggests.
best_label <- results$config[which.min(results$MAPE_test)]
cat("\nSelected:", best_label, "\n")

# ---- Step 2: refit selected config, forecast, evaluate ----------------
fit_tbats <- fits[[best_label]]
print(fit_tbats)

fc_tbats <- forecast(fit_tbats, h = H)
print(autoplot(fc_tbats) + autolayer(test, series = "Actual") +
        ggtitle(paste(best_label, "forecast")))

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_tbats, test))

# ---- Step 3: residual diagnostics -------------------------------------
cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_tbats)
print(Box.test(residuals(fit_tbats), lag = LAG_MAX, type = "Ljung-Box"))
cat("ACF-out-of-bounds:", acf_out_of_bounds(residuals(fit_tbats)), "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_tbats, test))
cat("Train/test MAPE gap:", round(g$mape_gap * 100, 1),
    "% | within 10%:", g$within_10pct, "\n")

# ---- Step 4: save in the shared summary shape -------------------------
summ <- model_summary("TBATS", fit_tbats, fc_tbats, test)
print(summ)
write.csv(summ, "output/summary_tbats.csv", row.names = FALSE)
saveRDS(fit_tbats, "output/fit_tbats.rds")
saveRDS(fc_tbats,  "output/fc_tbats.rds")
