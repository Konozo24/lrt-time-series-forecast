# 07_bats.R - Model family: Exponential Smoothing (BATS variant).
#
# BATS = Box-Cox transformation, ARMA errors, Trend, Seasonal components
# (De Livera, Hyndman & Snyder, 2011).
#
# Why BATS rather than TBATS here: TBATS replaces the seasonal states with
# trigonometric (Fourier) terms, a representation designed for COMPLEX
# seasonality - multiple seasonal periods, non-integer periods, or a very
# high frequency. This series has none of those: monthly ridership with a
# single integer period m = 12. With no complex seasonality to represent,
# the trigonometric machinery is not doing the job it was designed for,
# so BATS - which keeps the conventional seasonal-state representation -
# is the more direct specification for this data.
#
# The trade-off, stated honestly: BATS estimates a full set of 12 seasonal
# states where TBATS would use k <= 6 Fourier pairs, so BATS spends more
# parameters on a 79-observation training window. Whether that costs
# accuracy is exactly what the AIC comparison below and the holdout in
# 08_group_comparison.R are there to answer.
#
# Role in the comparison: an extended exponential-smoothing state-space
# model, tested against the conventional ETS in 06_ets.R. BATS adds to
# ETS a Box-Cox transform (to stabilise a seasonal swing that scales with
# the level), an ARMA correction on the residuals, and an optional damped
# trend; whether any of these is actually needed here is tested rather
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
  list(label = "BATS(box-cox, damped)",      bc = TRUE,  damped = TRUE),
  list(label = "BATS(box-cox, undamped)",    bc = TRUE,  damped = FALSE),
  list(label = "BATS(no box-cox, damped)",   bc = FALSE, damped = TRUE),
  list(label = "BATS(no box-cox, undamped)", bc = FALSE, damped = FALSE)
)

results <- data.frame()
fits <- list()

for (cfg in configs) {
  fit <- tryCatch(
    bats(train, use.box.cox = cfg$bc, use.trend = TRUE,
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
cat("\n--- BATS candidate configurations ---\n")
print(results, row.names = FALSE)
write.csv(results, tbl("bats_grid.csv"), row.names = FALSE)

# Selection is by AIC (results are already sorted by it). Holdout MAPE is
# printed alongside only as a cross-check, and is NOT the selection
# criterion - choosing a configuration by test-set accuracy would leak the
# test set into model selection.
best_label <- results$config[1]
cat("\nSelected (by AIC):", best_label, "\n")

# Whether AIC and holdout MAPE agree is worth reporting either way: if the
# AIC-best configuration is also the most accurate out of sample, no
# trade-off has to be argued; if they disagree, say so and keep the AIC
# choice, because switching to the MAPE-best one would leak the test set.
if (which.min(results$MAPE_test) != 1) {
  cat("NOTE: AIC and holdout MAPE disagree - AIC-best is",
      best_label, "but MAPE-best is",
      results$config[which.min(results$MAPE_test)],
      "\n  Keeping the AIC choice (selecting on test MAPE would leak the holdout).\n")
}

# ---- Step 2: refit selected config, forecast, evaluate ----------------
fit_bats <- fits[[best_label]]
print(fit_bats)

fc_bats <- forecast(fit_bats, h = H)
p_fc <- autoplot(fc_bats) + autolayer(test, series = "Actual") +
  ggtitle(paste(best_label, "forecast")) + scale_y_millions()
print(p_fc)
ggsave(fig("models", "bats_forecast.png"), p_fc, width = 8, height = 5, dpi = 150)

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_bats, test))

# ---- Step 3: residual diagnostics -------------------------------------
cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_bats)
dev.copy(png, filename = fig("models", "bats_residuals.png"), width = 900, height = 700, res = 130)
dev.off()
# Both lags: a spike at the seasonal lag stands out in a lag-12 test but
# gets diluted among mostly-zero terms in a longer one, so clear both.
# lb_test() supplies fitdf, which is 0 for BATS - matching how
# checkresiduals() treats it.
print(lb_test(fit_bats, 12))
print(lb_test(fit_bats, LAG_MAX))
cat("ACF-out-of-bounds: lag12 =", acf_out_of_bounds(residuals(fit_bats), lag.max = 12),
    "/ 12 | lag", LAG_MAX, "=", acf_out_of_bounds(residuals(fit_bats), lag.max = LAG_MAX),
    "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_bats, test))
cat("\nOverfitting checks (", g$direction, ")\n")
cat("  MASE gap  :", round(g$mase_gap * 100, 1), "% | <= 10%  :", g$within_10pct, "\n")
cat("  RMSE ratio:", round(g$rmse_ratio, 3), "  | <= 1.3x :", g$within_1_3x, "\n")
cat("  (MAPE gap", round(g$mape_gap * 100, 1),
    "% reported only - distorted by the level difference between the\n",
    "  MCO-era training window and the test window, so not part of the rule)\n")

# ---- Step 4: save in the shared summary shape -------------------------
summ <- model_summary("BATS", fit_bats, fc_bats, test)
print(summ)
write.csv(summ, tbl("summary_bats.csv"), row.names = FALSE)
