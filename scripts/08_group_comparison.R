# 08_group_comparison.R - combines all four models plus the SNAIVE
# baseline into one comparison on the same series, same split, same
# metrics. Run 04-07 first (this script refits them rather than reloading,
# so the comparison is self-contained and reproducible in one pass).
#
# SNAIVE is included as the benchmark, not as a fifth "real" model: it
# has no estimated parameters (it simply repeats the value from the same
# month one year earlier). It is the appropriate baseline here rather
# than a plain naive or mean forecast, because the series has confirmed
# seasonality (03_eda.R) - benchmarking against a method that ignores
# seasonality would be an artificially weak comparison. Any model that
# cannot beat SNAIVE is not earning the complexity it adds.
#
# Two comparisons are reported, and they answer different questions:
#   1. Against the baseline - is the modelling effort worth anything?
#   2. Within family (ARIMA vs SARIMA, ETS vs BATS) - does the seasonal
#      or extended specification earn its extra parameters?

source("scripts/00_setup.R")

y  <- load_series()
sp <- split_series(y)
train <- sp$train; test <- sp$test

cat("Series:", length(y), "months | train:", length(train),
    "| test:", length(test), "(h =", H, ")\n\n")

# ---- Refit all four models + baseline ---------------------------------
# Orders/specs below are the ones selected by the grid searches in 04-07.
# They are hard-coded here so this script runs standalone; if a grid
# search selects something different when re-run, update these to match.

fits <- list()

fits[["SNAIVE"]] <- snaive(train, h = H)

d_ns <- ndiffs(train, test = "kpss")
arima_grid <- expand.grid(p = 0:3, q = 0:3)
arima_res <- data.frame()
for (i in seq_len(nrow(arima_grid))) {
  f <- tryCatch(Arima(train, order = c(arima_grid$p[i], d_ns, arima_grid$q[i]),
                      include.drift = (d_ns > 0)), error = function(e) NULL)
  if (!is.null(f)) arima_res <- rbind(arima_res,
      data.frame(p = arima_grid$p[i], q = arima_grid$q[i], AICc = f$aicc))
}
arima_res <- arima_res[order(arima_res$AICc), ]
fits[["ARIMA"]] <- Arima(train, order = c(arima_res$p[1], d_ns, arima_res$q[1]),
                         include.drift = (d_ns > 0))

D_s <- nsdiffs(train)
d_s <- ndiffs(if (D_s > 0) diff(train, lag = 12) else train, test = "kpss")
sar_grid <- expand.grid(p = 0:2, q = 0:2, P = 0:1, Q = 0:1)
sar_res <- data.frame()
for (i in seq_len(nrow(sar_grid))) {
  f <- tryCatch(Arima(train, order = c(sar_grid$p[i], d_s, sar_grid$q[i]),
                      seasonal = list(order = c(sar_grid$P[i], D_s, sar_grid$Q[i]),
                                      period = 12)),
                error = function(e) NULL)
  if (!is.null(f)) sar_res <- rbind(sar_res, data.frame(
      p = sar_grid$p[i], q = sar_grid$q[i],
      P = sar_grid$P[i], Q = sar_grid$Q[i], AICc = f$aicc))
}
sar_res <- sar_res[order(sar_res$AICc), ]
fits[["SARIMA"]] <- Arima(train,
  order    = c(sar_res$p[1], d_s, sar_res$q[1]),
  seasonal = list(order = c(sar_res$P[1], D_s, sar_res$Q[1]), period = 12))

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

# Config matches the one selected in 07_bats.R by AIC (Box-Cox off,
# damped trend on) - see output/tables/bats_grid.csv. Do not change these
# flags without also updating 07_bats.R and 09_rolling_cv.R to match.
fits[["BATS"]] <- bats(train, use.box.cox = FALSE, use.trend = TRUE,
                       use.damped.trend = TRUE)

# ---- Build the comparison table ---------------------------------------
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
    direction    = g$direction,
    mape_gap_pct = round(g$mape_gap * 100, 1)))
}
rows <- rows[order(rows$MAPE_test), ]

cat("\n=== MODEL COMPARISON (ranked by test MAPE) ===\n")
print(rows[, c("model", "MAPE_test", "RMSE_test", "MAE_test", "MASE_test")],
      row.names = FALSE)

# ---- Overfitting checks, split out so both rules are legible ----------
# Rule 1: MASE gap <= 10%.  Rule 2: RMSE_test/RMSE_train <= 1.3.
# `direction` matters: on this series train error usually EXCEEDS test
# error, because the training window spans the MCO disruption and the
# recovery ramp while the test window is a flat mature period. That is
# the opposite of overfitting, so a large gap in the "train worse"
# direction is not a red flag.
cat("\n=== OVERFITTING CHECKS ===\n")
print(rows[, c("model", "MASE_train", "MASE_test", "mase_gap_pct", "within_10pct",
               "RMSE_train", "RMSE_test", "rmse_ratio", "within_1_3x", "direction")],
      row.names = FALSE)

cat("\n=== RESIDUAL DIAGNOSTICS (both lags, fitdf-corrected) ===\n")
print(rows[, c("model", "fitdf", "lb_pvalue_12", "lb_pvalue_16", "lb_pass",
               "n_lags_out_12", "n_lags_out_16")], row.names = FALSE)
cat("\nlb_pass = white noise at BOTH lags (p > 0.05). fitdf = p+q+P+Q, the\n",
    "parameters the ARIMA family estimates; ignoring it inflates p-values.\n")

# ---- Skill vs. the SNAIVE baseline ------------------------------------
# Skill score = 1 - (model error / baseline error). Positive means the
# model beats the benchmark; negative means it does not.
base_mape <- rows$MAPE_test[rows$model == "SNAIVE"]
rows$skill_vs_snaive <- round(1 - rows$MAPE_test / base_mape, 3)

cat("\n=== SKILL VS SNAIVE BASELINE ===\n")
print(rows[, c("model", "MAPE_test", "skill_vs_snaive")], row.names = FALSE)

write.csv(rows, tbl("model_comparison.csv"), row.names = FALSE)

# ---- Within-family comparison ----------------------------------------
cat("\n=== WITHIN-FAMILY COMPARISON ===\n")
cat("ARIMA family - does adding seasonal terms pay off?\n")
cat("  ARIMA  AICc:", round(fits[["ARIMA"]]$aicc, 2),
    " | SARIMA AICc:", round(fits[["SARIMA"]]$aicc, 2), "\n")
cat("ETS family - do BATS's Box-Cox / ARMA-error / damped-trend extensions\n")
cat("  earn their extra parameters over plain ETS?\n")
cat("  ETS test MAPE:",   rows$MAPE_test[rows$model == "ETS"],
    " | BATS test MAPE:", rows$MAPE_test[rows$model == "BATS"], "\n")

# ---- Combined forecast plot ------------------------------------------
fc_all <- lapply(names(fits), function(nm)
  if (nm == "SNAIVE") fits[[nm]] else forecast(fits[[nm]], h = H))
names(fc_all) <- names(fits)

p <- autoplot(window(y, start = c(2022, 1))) +
  autolayer(fc_all[["ARIMA"]]$mean,  series = "ARIMA") +
  autolayer(fc_all[["SARIMA"]]$mean, series = "SARIMA") +
  autolayer(fc_all[["ETS"]]$mean,    series = "ETS") +
  autolayer(fc_all[["BATS"]]$mean,   series = "BATS") +
  autolayer(fc_all[["SNAIVE"]]$mean, series = "SNAIVE") +
  autolayer(test, series = "Actual", size = 1.1) +
  ggtitle("All models vs. actual (12-month holdout)") +
  ylab("Monthly ridership") + guides(colour = guide_legend(title = "Model")) +
  scale_y_millions()
print(p)
ggsave(fig("comparison", "model_comparison_plot.png"), p, width = 10, height = 6, dpi = 150)

# ---- Bar chart: MAPE by model, benchmark called out ------------------
p_bar <- ggplot(rows, aes(x = reorder(model, MAPE_test), y = MAPE_test,
                          fill = model == "SNAIVE")) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(MAPE_test, "%")), hjust = -0.15, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "grey60", "FALSE" = "steelblue")) +
  labs(title = "Holdout MAPE by model (SNAIVE benchmark in grey)",
       x = NULL, y = "MAPE (%)") +
  theme_minimal()
print(p_bar)
ggsave(fig("comparison", "model_comparison_bar.png"), p_bar, width = 8, height = 5, dpi = 150)
