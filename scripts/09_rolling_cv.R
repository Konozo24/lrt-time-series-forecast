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
# RUNTIME WARNING: this refits every model at every origin, and BATS is
# slow because it runs its own internal component search on each fit.
# Expect several minutes on free Posit Cloud. Reduce the number of folds
# (widen the `by` argument) if it becomes impractical.

source("scripts/00_setup.R")

y <- load_series()
n <- length(y)
h <- H                      # 12-month horizon per fold

origins <- seq(n - h - 16, n - h, by = 4)   # 5 expanding-window origins
cat("Folds:", length(origins), "| training sizes:", paste(origins, collapse = ", "), "\n\n")

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

  # ARIMA (non-seasonal), order re-selected per fold by AICc
  d1 <- ndiffs(tr, test = "kpss")
  gA <- expand.grid(p = 0:2, q = 0:2); rA <- data.frame()
  for (i in seq_len(nrow(gA))) {
    f <- tryCatch(Arima(tr, order = c(gA$p[i], d1, gA$q[i]),
                        include.drift = (d1 > 0)), error = function(e) NULL)
    if (!is.null(f)) rA <- rbind(rA, data.frame(p = gA$p[i], q = gA$q[i], AICc = f$aicc))
  }
  rA <- rA[order(rA$AICc), ]
  m_ar <- metrics_of(forecast(Arima(tr, order = c(rA$p[1], d1, rA$q[1]),
                                    include.drift = (d1 > 0)), h = h)$mean)

  # SARIMA, orders re-selected per fold
  D1 <- nsdiffs(tr)
  d2 <- ndiffs(if (D1 > 0) diff(tr, lag = 12) else tr, test = "kpss")
  gS <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1); rS <- data.frame()
  for (i in seq_len(nrow(gS))) {
    f <- tryCatch(Arima(tr, order = c(gS$p[i], d2, gS$q[i]),
                        seasonal = list(order = c(gS$P[i], D1, gS$Q[i]), period = 12)),
                  error = function(e) NULL)
    if (!is.null(f)) rS <- rbind(rS, data.frame(p = gS$p[i], q = gS$q[i],
                                                P = gS$P[i], Q = gS$Q[i], AICc = f$aicc))
  }
  rS <- rS[order(rS$AICc), ]
  m_sa <- metrics_of(forecast(Arima(tr, order = c(rS$p[1], d2, rS$q[1]),
                                    seasonal = list(order = c(rS$P[1], D1, rS$Q[1]),
                                                    period = 12)), h = h)$mean)

  # ETS, spec re-selected per fold by AICc
  specs <- list(c("ANN", NA), c("AAN", "TRUE"), c("ANA", NA),
                c("AAA", "TRUE"), c("MNM", NA), c("MAM", "TRUE"))
  rE <- data.frame(); eF <- list()
  for (s in specs) {
    dmp <- if (is.na(s[2])) NULL else as.logical(s[2])
    f <- tryCatch(ets(tr, model = s[1], damped = dmp), error = function(e) NULL)
    if (!is.null(f)) {
      k <- paste0(s[1], ifelse(is.na(s[2]), "", s[2])); eF[[k]] <- f
      rE <- rbind(rE, data.frame(key = k, AICc = f$aicc))
    }
  }
  rE <- rE[order(rE$AICc), ]
  m_ets <- metrics_of(forecast(eF[[rE$key[1]]], h = h)$mean)

  # BATS - config as selected in 07_bats.R (Box-Cox off, damped trend on)
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
