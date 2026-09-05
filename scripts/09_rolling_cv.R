# 09_rolling_cv.R - rolling-origin cross-validation for all four models.
#
# Why this exists: the comparison in 08_group_comparison.R rests on a
# single 12-month holdout. One holdout is one draw - if that particular
# window happens to be atypical (an unusually flat or volatile year), the
# ranking it produces may not hold. Rolling-origin CV refits each model
# at several successive origins and averages accuracy across them, so the
# ranking reflects consistent performance rather than one lucky window.
#
# Design: expanding training window, 12-month forecast horizon per fold
# (one full seasonal cycle, matching the main holdout), 5 folds.
#
# Configs: SNAIVE/SARIMA/ETS/BATS are refit at each origin using the SAME
# hard-coded config as 08_group_comparison.R (each model's own grid-search
# winner on the full 79-month training set), rather than re-searching a
# grid at every fold. This is a deliberate trade-off, not an oversight:
# refitting the already-chosen structure at each origin still measures
# genuine out-of-sample accuracy (the parameters ARE re-estimated per
# fold, just not the order/spec), and it keeps CV fast. What this does
# NOT measure is whether a shorter training window would have picked a
# different order/spec - a fold-specific model this script does not
# search for.
#
# ARIMA is the exception: it's a STLARIMA (see 00_setup.R's
# stlarima_fit()/stlarima_forecast() and 04_arima.R), and a bare order
# tuple isn't a complete spec for it - it has to be paired with each
# fold's own STL decomposition and xreg. So its order IS re-searched per
# fold here (reduced grid, p/q = 0:2), matching what 04_arima.R's own
# rolling-CV section already does for the same reason. Do not change
# these without also updating 08_group_comparison.R and 04-07 together.

source("scripts/00_setup.R")

y <- load_series()
n <- length(y)
h <- H                      # 12-month horizon per fold

origins <- seq(n - h - 16, n - h, by = 4)   # 5 expanding-window origins
cat("Folds:", length(origins), "| training sizes:", paste(origins, collapse = ", "), "\n\n")

xreg_all <- intervention_xreg(y)   # STLARIMA (ARIMA row) intervention pulses

cv <- data.frame()

for (ts_size in origins) {
  tr <- head(y, ts_size)
  te <- window(y, start = time(y)[ts_size + 1], end = time(y)[ts_size + h])

  # Three metrics per fold, not just MAPE. RMSE is needed for the
  # RMSE_test/RMSE_train <= 1.3 rule judged against the CV mean (the
  # single holdout is one draw; the CV mean is the more reliable
  # denominator), and MASE is the scale-free accuracy measure. MASE is
  # scaled by THIS fold's in-sample seasonal-naive MAE, matching
  # forecast::accuracy()'s definition for seasonal data.
  metrics_of <- function(fc_mean) {
    e     <- as.numeric(te) - as.numeric(fc_mean)
    scale <- mean(abs(diff(as.numeric(tr), lag = 12)))
    c(MAPE = mean(abs(e / as.numeric(te))) * 100,
      RMSE = sqrt(mean(e^2)),
      MASE = mean(abs(e)) / scale)
  }
  # SNAIVE baseline
  m_sn <- metrics_of(snaive(tr, h = h)$mean)

  # STLARIMA - order re-searched on THIS fold's own STL decomposition
  # (reduced grid, p/q = 0:2), not fixed to the full-series winner - see
  # header. xreg pulses come from this fold's own train/test slice.
  xreg_tr <- xreg_all[1:ts_size, , drop = FALSE]
  xreg_te <- xreg_all[(ts_size + 1):(ts_size + h), , drop = FALSE]
  sf_ar   <- stlarima_fit(tr, xreg_tr, p_range = 0:2, q_range = 0:2)
  m_ar    <- metrics_of(stlarima_forecast(sf_ar, tr, h, xreg_te)$mean)

  # SARIMA(0,1,2)(1,0,1)[12] - config fixed to the full-series winner.
  m_sa <- metrics_of(forecast(Arima(tr, order = c(0, 1, 2),
                                    seasonal = list(order = c(1, 0, 1), period = 12)),
                              h = h)$mean)

  # ETS(A,N,A) - config fixed to the full-series winner.
  m_ets <- metrics_of(forecast(ets(tr, model = "ANA", damped = NULL), h = h)$mean)

  # BATS(no box-cox, damped trend) - config fixed to the full-series winner.
  m_tb <- metrics_of(forecast(bats(tr, use.box.cox = FALSE, use.trend = TRUE,
                                    use.damped.trend = TRUE), h = h)$mean)

  # Long format: one row per model per fold, three metrics each.
  fold <- rbind(m_sn, m_ar, m_sa, m_ets, m_tb)
  cv <- rbind(cv, data.frame(train_size = ts_size,
                             model = c("SNAIVE", "ARIMA", "SARIMA", "ETS", "BATS"),
                             MAPE = fold[, "MAPE"], RMSE = fold[, "RMSE"],
                             MASE = fold[, "MASE"], row.names = NULL))
  cat("origin", ts_size, "done\n")
}

models <- c("SNAIVE", "ARIMA", "SARIMA", "ETS", "BATS")

cat("\n=== Per-fold MAPE (%) ===\n")
# one value per (fold, model) cell, so mean() just unwraps it into a
# numeric matrix - tapply with identity can come back as a list.
print(round(with(cv, tapply(MAPE, list(train_size, model), mean))[, models], 2))

summary_cv <- do.call(rbind, lapply(models, function(m) {
  d <- cv[cv$model == m, ]
  data.frame(model     = m,
             mean_MAPE = round(mean(d$MAPE), 2),
             sd_MAPE   = round(sd(d$MAPE), 2),
             min_MAPE  = round(min(d$MAPE), 2),
             max_MAPE  = round(max(d$MAPE), 2),
             mean_RMSE = round(mean(d$RMSE), 0),
             mean_MASE = round(mean(d$MASE), 3),
             n_folds   = nrow(d))
}))
summary_cv <- summary_cv[order(summary_cv$mean_MAPE), ]

cat("\n=== Rolling-CV summary (ranked by mean MAPE) ===\n")
print(summary_cv, row.names = FALSE)
cat("\nsd_MAPE is the stability measure: a model with a low mean but high sd\n",
    "performs well on average while being unreliable fold to fold.\n")

write.csv(cv, tbl("rolling_cv_folds.csv"), row.names = FALSE)
write.csv(summary_cv, tbl("rolling_cv_summary.csv"), row.names = FALSE)

# ---- Overfitting checks judged against the CV mean --------------------
# 08_group_comparison.R applies both rules against the SINGLE holdout.
# That holdout is one draw, so the same rules are re-applied here with the
# CV mean as the test-side term - the more reliable version, and the one
# to quote in the report. Training-side terms come from the full-training
# fits in model_comparison.csv, so run 08 before this script.
main_path <- tbl("model_comparison.csv")
if (file.exists(main_path)) {
  main <- read.csv(main_path, stringsAsFactors = FALSE)
  if (all(c("RMSE_train", "MASE_train") %in% names(main))) {
    ov <- merge(summary_cv[, c("model", "mean_RMSE", "mean_MASE")],
                main[, c("model", "RMSE_train", "MASE_train",
                         "RMSE_test", "MASE_test")], by = "model")
    ov$rmse_ratio_cv      <- round(ov$mean_RMSE / ov$RMSE_train, 3)
    ov$within_1_3x_cv     <- ov$rmse_ratio_cv <= 1.3
    ov$gap_pct_cv         <- round(abs(ov$mean_MASE - ov$MASE_train) /
                                     ov$mean_MASE * 100, 1)
    ov$within_10pct_cv    <- ov$gap_pct_cv <= 10
    ov$rmse_ratio_holdout <- round(ov$RMSE_test / ov$RMSE_train, 3)
    ov <- ov[order(ov$rmse_ratio_cv), ]

    cat("\n=== OVERFITTING CHECKS vs CV MEAN (authoritative) ===\n")
    print(ov[, c("model", "RMSE_train", "mean_RMSE", "rmse_ratio_cv",
                 "within_1_3x_cv", "MASE_train", "mean_MASE", "gap_pct_cv",
                 "within_10pct_cv", "rmse_ratio_holdout")], row.names = FALSE)
    cat("\nrmse_ratio_holdout is shown for reference only - the CV column is the\n",
        "one to report, since it averages over five origins instead of one.\n")
    write.csv(ov, tbl("overfit_checks_cv.csv"), row.names = FALSE)
  } else {
    cat("\nNote: model_comparison.csv has no RMSE_train/MASE_train columns.\n",
        "Re-run 08_group_comparison.R to regenerate it, then re-run this script.\n")
  }
} else {
  cat("\nNote:", main_path, "not found - run 08_group_comparison.R first to get\n",
      "the CV-based overfitting checks.\n")
}

# ---- Stability plot: mean MAPE vs. sd MAPE -----------------------------
# The point of rolling CV is this trade-off, not just the mean ranking:
# a model can look good on average (low x) but be unreliable fold-to-fold
# (high y). The best model is bottom-left; SNAIVE is expected top-right.
p_stab <- ggplot(summary_cv, aes(x = mean_MAPE, y = sd_MAPE, label = model)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -1, size = 4) +
  labs(title = "Rolling-CV: accuracy vs. stability (bottom-left is best)",
       x = "Mean MAPE across folds (%)", y = "SD of MAPE across folds (stability)") +
  theme_minimal()
print(p_stab)
ggsave(fig("comparison", "rolling_cv_stability.png"), p_stab, width = 7, height = 5, dpi = 150)

# ---- Per-fold line plot: does any model's rank flip across folds? -----
# `cv` is already long (one row per model per fold), so it plots directly.
p_folds <- ggplot(cv, aes(x = train_size, y = MAPE, color = model)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  labs(title = "MAPE per rolling-CV fold, by model",
       x = "Training window size (months)", y = "MAPE (%)") +
  theme_minimal()
print(p_folds)
ggsave(fig("comparison", "rolling_cv_per_fold.png"), p_folds, width = 8, height = 5, dpi = 150)
