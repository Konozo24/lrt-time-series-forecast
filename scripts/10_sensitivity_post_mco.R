# 10_sensitivity_post_mco.R - does the MCO reconstruction drive the result?
#
# THE QUESTION THIS ANSWERS: 22 of the 79 training months used in scripts
# 04-09 are not observations - they were reconstructed in 02_mco_resolution.R
# (bridged trend + STL seasonal). That is 28% of the training data. So a
# fair challenge to the whole project is: did BATS win because it is the
# better model, or because it happened to fit the constructed segment well?
#
# THE TEST: refit the same models on the post-MCO subsample only - 2022-01
# onward - where every single observation is real and nothing is imputed.
# Then compare the ranking against the main analysis.
#   - Same winner  -> the reconstruction is not driving the conclusion, and
#                     that can be stated in the report with evidence.
#   - Different winner -> the reconstruction WAS driving it, which is far
#                     better to discover here than to have it found later.
#
# PRECEDENT: this is the same device used by Hightower et al. (2024,
# Journal of Public Transportation 26, 100097), who estimated their models
# over three separate periods (pre-COVID, full series, post-COVID) rather
# than adjusting the data, specifically to isolate the pandemic's effect on
# forecast performance.
#
# READ THE OUTPUT WITH THIS CAVEAT: only 55 months are available from
# 2022-01, so training is 43 months = 3.6 seasonal cycles, against 6.6 in
# the main analysis. That is thin for a seasonal model - SARIMA in
# particular - so absolute errors here are expected to be worse and less
# stable than the main results. This script is a ROBUSTNESS CHECK on the
# RANKING, not a replacement analysis, and it is not a fair head-to-head
# on accuracy.

source("scripts/00_setup.R")

# Deliberately load the RAW series, not the resolved one, so there is no
# question of imputed values leaking in. (The two are identical from
# 2022-01 onward - only the 22 MCO months ever differed - but loading the
# raw file makes the "zero reconstructed observations" claim unambiguous.)
y_raw <- readRDS("data/ampang_monthly_full.rds")
y <- window(y_raw, start = c(2022, 1))

sp <- split_series(y)          # same h = 12 holdout as the main analysis
train <- sp$train; test <- sp$test

cat("Post-MCO subsample:", length(y), "months (",
    format(zoo::as.yearmon(time(y)[1])), "to",
    format(zoo::as.yearmon(time(y)[length(y)])), ")\n")
cat("Train:", length(train), "months =",
    round(length(train) / 12, 2), "seasonal cycles | test:", length(test), "\n")
cat("Reconstructed observations in this subsample: 0\n\n")

fits <- list()

# ---- Baseline ---------------------------------------------------------
fits[["SNAIVE"]] <- snaive(train, h = H)

# ---- ARIMA: same grid and same selection rule as 04_arima.R -----------
d_ns <- ndiffs(train, test = "kpss")
arima_res <- data.frame()
for (i in seq_len(nrow(g <- expand.grid(p = 0:3, q = 0:3)))) {
  f <- tryCatch(Arima(train, order = c(g$p[i], d_ns, g$q[i]),
                      include.drift = (d_ns > 0)), error = function(e) NULL)
  if (!is.null(f)) arima_res <- rbind(arima_res,
      data.frame(p = g$p[i], q = g$q[i], AICc = f$aicc))
}
arima_res <- arima_res[order(arima_res$AICc), ]
fits[["ARIMA"]] <- Arima(train, order = c(arima_res$p[1], d_ns, arima_res$q[1]),
                         include.drift = (d_ns > 0))

# ---- SARIMA: same grid as 05_sarima.R ---------------------------------
# nsdiffs() needs enough cycles to estimate seasonal strength; on 43 months
# it can fail outright, so fall back to D = 0 rather than aborting the run.
D_s <- tryCatch(nsdiffs(train), error = function(e) 0)
d_s <- ndiffs(if (D_s > 0) diff(train, lag = 12) else train, test = "kpss")
sar_res <- data.frame()
for (i in seq_len(nrow(g <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)))) {
  f <- tryCatch(Arima(train, order = c(g$p[i], d_s, g$q[i]),
                      seasonal = list(order = c(g$P[i], D_s, g$Q[i]), period = 12)),
                error = function(e) NULL)
  if (!is.null(f)) sar_res <- rbind(sar_res, data.frame(
      p = g$p[i], q = g$q[i], P = g$P[i], Q = g$Q[i], AICc = f$aicc))
}
sar_res <- sar_res[order(sar_res$AICc), ]
fits[["SARIMA"]] <- Arima(train,
  order    = c(sar_res$p[1], d_s, sar_res$q[1]),
  seasonal = list(order = c(sar_res$P[1], D_s, sar_res$Q[1]), period = 12))

# ---- ETS: same candidate list as 06_ets.R -----------------------------
ets_specs <- list(c("ANN", NA), c("AAN", "FALSE"), c("AAN", "TRUE"),
                  c("ANA", NA), c("AAA", "FALSE"), c("AAA", "TRUE"),
                  c("MNM", NA), c("MAM", "FALSE"), c("MAM", "TRUE"))
ets_res <- data.frame(); ets_fits <- list()
for (s in ets_specs) {
  dmp <- if (is.na(s[2])) NULL else as.logical(s[2])
  f <- tryCatch(ets(train, model = s[1], damped = dmp), error = function(e) NULL)
  if (!is.null(f)) {
    key <- paste0(s[1], "_", ifelse(is.na(s[2]), "NA", s[2]))
    ets_fits[[key]] <- f
    ets_res <- rbind(ets_res, data.frame(key = key, AICc = f$aicc))
  }
}
ets_res <- ets_res[order(ets_res$AICc), ]
fits[["ETS"]] <- ets_fits[[ets_res$key[1]]]

# ---- BATS: components RE-SELECTED here, not inherited ------------------
# The main analysis picked its config by AIC on the full 79-month training
# set. Reusing that choice here would import a decision made with the help
# of the reconstructed months - exactly what this script is meant to avoid.
bats_cfgs <- list(list(bc = TRUE,  damped = TRUE),  list(bc = TRUE,  damped = FALSE),
                  list(bc = FALSE, damped = TRUE),  list(bc = FALSE, damped = FALSE))
bats_res <- data.frame(); bats_fits <- list()
for (k in seq_along(bats_cfgs)) {
  cfg <- bats_cfgs[[k]]
  f <- tryCatch(bats(train, use.box.cox = cfg$bc, use.trend = TRUE,
                     use.damped.trend = cfg$damped), error = function(e) NULL)
  if (!is.null(f)) {
    lab <- sprintf("BATS(%sbox-cox, %s)", if (cfg$bc) "" else "no ",
                   if (cfg$damped) "damped" else "undamped")
    bats_fits[[lab]] <- f
    bats_res <- rbind(bats_res, data.frame(config = lab, AIC = f$AIC))
  }
}
bats_res <- bats_res[order(bats_res$AIC), ]
cat("--- BATS configs re-selected on the post-MCO subsample ---\n")
print(bats_res, row.names = FALSE)
cat("Selected (by AIC):", bats_res$config[1], "\n\n")
fits[["BATS"]] <- bats_fits[[bats_res$config[1]]]

# ---- Comparison table, identical shape to 08_group_comparison.R --------
rows <- data.frame()
for (nm in names(fits)) {
  fit <- fits[[nm]]
  fc  <- if (nm == "SNAIVE") fit else forecast(fit, h = H)
  a   <- accuracy(fc, test)
  g   <- gap_check(a)
  r   <- residuals(fit)
  p12 <- lb_test(fit, 12)$p.value
  p16 <- lb_test(fit, LAG_MAX)$p.value
  rows <- rbind(rows, data.frame(
    model        = nm,
    MAPE_test    = round(a["Test set", "MAPE"], 3),
    RMSE_test    = round(a["Test set", "RMSE"], 0),
    MAE_test     = round(a["Test set", "MAE"], 0),
    MASE_test    = round(a["Test set", "MASE"], 3),
    RMSE_train   = round(g$rmse_train, 0),
    MASE_train   = round(g$mase_train, 3),
    fitdf        = model_fitdf(fit),
    lb_pvalue_12 = round(p12, 4),
    lb_pvalue_16 = round(p16, 4),
    lb_pass      = (p12 > 0.05) & (p16 > 0.05),
    n_lags_out_12 = acf_out_of_bounds(r, lag.max = 12),
    n_lags_out_16 = acf_out_of_bounds(r, lag.max = LAG_MAX),
    mase_gap_pct = round(g$mase_gap * 100, 1),
    within_10pct = g$within_10pct,
    rmse_ratio   = round(g$rmse_ratio, 3),
    within_1_3x  = g$within_1_3x,
    direction    = g$direction))
}
rows <- rows[order(rows$MAPE_test), ]

cat("=== POST-MCO SUBSAMPLE: model comparison (ranked by test MAPE) ===\n")
print(rows, row.names = FALSE)
write.csv(rows, tbl("sensitivity_post_mco.csv"), row.names = FALSE)

# ---- The actual verdict: did the ranking survive? ---------------------
# Compare against the main analysis rather than making the reader do it.
main_path <- tbl("model_comparison.csv")
if (file.exists(main_path)) {
  main <- read.csv(main_path, stringsAsFactors = FALSE)
  cmp <- merge(
    data.frame(model = main$model,  rank_full = seq_len(nrow(main)),
               MAPE_full = main$MAPE_test),
    data.frame(model = rows$model,  rank_post = seq_len(nrow(rows)),
               MAPE_post = rows$MAPE_test),
    by = "model")
  cmp <- cmp[order(cmp$rank_full), ]

  cat("\n=== RANKING: full series (91 months, 22 reconstructed)",
      "vs post-MCO (55 months, 0 reconstructed) ===\n")
  print(cmp, row.names = FALSE)

  win_full <- main$model[1]; win_post <- rows$model[1]
  cat("\nBest on full series:", win_full, "| best on post-MCO subsample:", win_post, "\n")
  if (identical(win_full, win_post)) {
    cat("VERDICT: same winner. The ranking does not depend on the reconstructed\n",
        " months - the conclusion holds on real observations alone.\n")
  } else {
    cat("VERDICT: the winner CHANGES.", win_full, "wins on the full series but",
        win_post, "wins\n on unimputed data. Report this openly: it means the\n",
        " reconstruction influences model selection, and the full-series ranking\n",
        " should be presented with that caveat rather than as a clean result.\n")
  }
  cat("\nSpearman rank correlation (full vs post-MCO):",
      round(cor(cmp$rank_full, cmp$rank_post, method = "spearman"), 3),
      "\n(1.00 = identical ordering; a high value with a changed winner still\n",
      " indicates broadly consistent model behaviour.)\n")
} else {
  cat("\nNote:", main_path, "not found - run 08_group_comparison.R first to get\n",
      "the side-by-side ranking comparison.\n")
}

# ---- Plot: both periods' forecasts against actuals --------------------
fc_all <- lapply(names(fits), function(nm)
  if (nm == "SNAIVE") fits[[nm]] else forecast(fits[[nm]], h = H))
names(fc_all) <- names(fits)

p <- autoplot(y) +
  autolayer(fc_all[["ARIMA"]]$mean,  series = "ARIMA") +
  autolayer(fc_all[["SARIMA"]]$mean, series = "SARIMA") +
  autolayer(fc_all[["ETS"]]$mean,    series = "ETS") +
  autolayer(fc_all[["BATS"]]$mean,   series = "BATS") +
  autolayer(fc_all[["SNAIVE"]]$mean, series = "SNAIVE") +
  autolayer(test, series = "Actual", linewidth = 1.1) +
  ggtitle("Sensitivity check: post-MCO subsample only (no reconstructed months)") +
  ylab("Monthly ridership") + guides(colour = guide_legend(title = "Model")) +
  scale_y_millions()
print(p)
ggsave(fig("comparison", "sensitivity_post_mco.png"), p, width = 10, height = 6, dpi = 150)

cat("\nWrote:", tbl("sensitivity_post_mco.csv"), "and",
    fig("comparison", "sensitivity_post_mco.png"), "\n")
