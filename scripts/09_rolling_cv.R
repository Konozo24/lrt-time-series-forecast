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
# RUNTIME WARNING: this refits every model at every origin, and TBATS is
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

  mape_of <- function(fc_mean) mean(abs((te - fc_mean) / te)) * 100

  # SNAIVE baseline
  m_sn <- mape_of(snaive(tr, h = h)$mean)

  # ARIMA (non-seasonal), order re-selected per fold by AICc
  d1 <- ndiffs(tr, test = "kpss")
  gA <- expand.grid(p = 0:2, q = 0:2); rA <- data.frame()
  for (i in seq_len(nrow(gA))) {
    f <- tryCatch(Arima(tr, order = c(gA$p[i], d1, gA$q[i]),
                        include.drift = (d1 > 0)), error = function(e) NULL)
    if (!is.null(f)) rA <- rbind(rA, data.frame(p = gA$p[i], q = gA$q[i], AICc = f$aicc))
  }
  rA <- rA[order(rA$AICc), ]
  m_ar <- mape_of(forecast(Arima(tr, order = c(rA$p[1], d1, rA$q[1]),
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
  m_sa <- mape_of(forecast(Arima(tr, order = c(rS$p[1], d2, rS$q[1]),
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
  m_ets <- mape_of(forecast(eF[[rE$key[1]]], h = h)$mean)

  # TBATS - config as selected in 07_tbats.R (Box-Cox off, damped trend on)
  m_tb <- mape_of(forecast(tbats(tr, use.box.cox = FALSE, use.trend = TRUE,
                                 use.damped.trend = TRUE), h = h)$mean)

  cv <- rbind(cv, data.frame(train_size = ts_size, SNAIVE = m_sn, ARIMA = m_ar,
                             SARIMA = m_sa, ETS = m_ets, TBATS = m_tb))
  cat("origin", ts_size, "done\n")
}

cat("\n=== Per-fold MAPE (%) ===\n")
print(round(cv, 2), row.names = FALSE)

summary_cv <- data.frame(
  model    = c("SNAIVE", "ARIMA", "SARIMA", "ETS", "TBATS"),
  mean_MAPE = round(sapply(cv[, -1], mean), 2),
  sd_MAPE   = round(sapply(cv[, -1], sd), 2),
  min_MAPE  = round(sapply(cv[, -1], min), 2),
  max_MAPE  = round(sapply(cv[, -1], max), 2))
summary_cv <- summary_cv[order(summary_cv$mean_MAPE), ]

cat("\n=== Rolling-CV summary (ranked by mean MAPE) ===\n")
print(summary_cv, row.names = FALSE)
cat("\nsd_MAPE is the stability measure: a model with a low mean but high sd\n",
    "performs well on average while being unreliable fold to fold.\n")

write.csv(cv, "output/rolling_cv_folds.csv", row.names = FALSE)
write.csv(summary_cv, "output/rolling_cv_summary.csv", row.names = FALSE)
