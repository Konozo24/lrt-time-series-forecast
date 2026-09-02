# 06_ets.R - Model family: Exponential Smoothing (ETS).
#
# Role in the comparison: the conventional exponential-smoothing member.
# Where ARIMA/SARIMA handle seasonality by differencing, ETS models it as
# explicit seasonal states updated by smoothing - a structurally different
# mechanism, which is why both families are worth comparing on the same
# series rather than assuming one representation suits the data.
#
# NOT using the automatic ets() model search: each candidate specification
# is fitted explicitly by name and ranked by AICc, so the comparison
# between additive and multiplicative seasonality (and between damped and
# undamped trend) is visible in the output rather than hidden inside an
# automatic selection step.
#
# ETS notation is ETS(Error, Trend, Seasonal):
#   Error    A = additive, M = multiplicative
#   Trend    N = none, A = additive, Ad = additive damped
#   Seasonal N = none, A = additive, M = multiplicative
# The additive-vs-multiplicative seasonal comparison is the substantive
# one here: multiplicative means the seasonal swing scales with the level
# of ridership rather than staying a fixed number of trips.

source("scripts/00_setup.R")

y  <- load_series()
sp <- split_series(y)
train <- sp$train; test <- sp$test

# ---- Step 1: explicit candidate specifications ------------------------
# Combinations are restricted to numerically stable ones: forecast::ets()
# treats additive error with multiplicative seasonality (A,*,M) as
# unstable, so those are fitted with multiplicative error instead.
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

results <- data.frame()
fits <- list()

for (cand in candidates) {
  fit <- tryCatch(
    ets(train, model = cand$model, damped = cand$damped),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    fits[[cand$label]] <- fit
    results <- rbind(results, data.frame(spec = cand$label,
                                         AICc = fit$aicc, BIC = fit$bic))
  }
}

results <- results[order(results$AICc), ]
cat("\n--- ETS candidate specifications (ranked by AICc) ---\n")
print(results, row.names = FALSE)
write.csv(results, tbl("ets_grid.csv"), row.names = FALSE)

best_label <- results$spec[1]
cat("\nSelected:", best_label, " AICc =", round(results$AICc[1], 2), "\n")

# ---- Step 2: refit selected spec, forecast, evaluate ------------------
fit_ets <- fits[[best_label]]
print(summary(fit_ets))

fc_ets <- forecast(fit_ets, h = H)
p_fc <- autoplot(fc_ets) + autolayer(test, series = "Actual") +
  ggtitle(paste(best_label, "forecast"))
print(p_fc)
ggsave(fig("models", "ets_forecast.png"), p_fc, width = 8, height = 5, dpi = 150)

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_ets, test))

# ---- Step 3: residual diagnostics -------------------------------------
cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_ets)
dev.copy(png, filename = fig("models", "ets_residuals.png"), width = 900, height = 700, res = 130)
dev.off()
# Both lags: a spike at the seasonal lag stands out in a lag-12 test but
# gets diluted among mostly-zero terms in a longer one, so clear both.
# lb_test() supplies fitdf, which is 0 for ETS - matching how
# checkresiduals() treats exponential-smoothing fits.
print(lb_test(fit_ets, 12))
print(lb_test(fit_ets, LAG_MAX))
cat("ACF-out-of-bounds: lag12 =", acf_out_of_bounds(residuals(fit_ets), lag.max = 12),
    "/ 12 | lag", LAG_MAX, "=", acf_out_of_bounds(residuals(fit_ets), lag.max = LAG_MAX),
    "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_ets, test))
cat("\nOverfitting checks (", g$direction, ")\n")
cat("  MASE gap  :", round(g$mase_gap * 100, 1), "% | <= 10%  :", g$within_10pct, "\n")
cat("  RMSE ratio:", round(g$rmse_ratio, 3), "  | <= 1.3x :", g$within_1_3x, "\n")
cat("  (MAPE gap", round(g$mape_gap * 100, 1),
    "% reported only - distorted by the level difference between the\n",
    "  MCO-era training window and the test window, so not part of the rule)\n")

# ---- Step 4: save in the shared summary shape -------------------------
summ <- model_summary(best_label, fit_ets, fc_ets, test)
print(summ)
write.csv(summ, tbl("summary_ets.csv"), row.names = FALSE)
