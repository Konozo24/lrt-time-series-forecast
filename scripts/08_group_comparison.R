# 08_group_comparison.R - combines all four models plus the SNAIVE
# baseline into one comparison on the same series, same split, same
# metrics. Refits all four itself rather than reloading anything - it
# does not read 04_arima.R/05_sarima.R/06_ets.R/07_bats.R's output (those
# four are standalone scripts, not part of this pipeline; see README) -
# so this script only needs 02_mco_resolution.R to have run first.
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
# SARIMA/ETS/BATS below are hard-coded to the winner of each model's own
# grid search (run in 05_sarima.R/06_ets.R/07_bats.R against this same
# 79-month training set) rather than re-searched here. This keeps every
# run fast and the comparison table stable and reproducible. Trade-off,
# stated honestly: if the source data changes (e.g. a new month is
# appended), these are NOT re-verified as still AICc-best - the
# corresponding standalone script must be re-run and this block updated
# to match. Do not change these without also updating 07_bats.R,
# 05_sarima.R, 06_ets.R, and 09_rolling_cv.R together.
#
# ARIMA is the one exception, and deliberately not hard-coded: 04_arima.R
# fits a STLARIMA (STL-decomposed ARIMA + intervention xreg, credited to
# Justin, github.com/JusunF/ARIMA/combine.R) rather than a plain ARIMA(p,d,q)
# on the raw series. A bare order tuple isn't a complete spec for that
# model - it has to be paired with the STL decomposition and the xreg it
# was selected against - so the shared stlarima_fit()/stlarima_forecast()
# helpers in 00_setup.R re-derive it here instead, on the same training
# set and with the same AICc-only, no-test-leakage selection 04_arima.R
# uses. See 00_setup.R for the helpers and 04_arima.R for the rationale.

fits <- list()

fits[["SNAIVE"]] <- snaive(train, h = H)

# STLARIMA - STL-decompose train, fit ARIMA(p,d,q)+xreg to the seasonally
# adjusted series (d by KPSS, p/q by AICc over {0:3, 0:3}), pulse xreg at
# the two intervention dates in 04_arima.R. fitdf = p+q (xreg coefficients
# don't count, matching model_fitdf()'s convention).
xreg_all   <- intervention_xreg(y)
xreg_train <- head(xreg_all, length(train))
xreg_test  <- tail(xreg_all, H)
sf_arima   <- stlarima_fit(train, xreg_train)
fits[["ARIMA"]] <- sf_arima$fit
fc_arima_precomputed <- stlarima_forecast(sf_arima, train, H, xreg_test)
cat("STLARIMA order:", paste(sf_arima$order, collapse = ","), "\n")

# SARIMA(0,1,2)(1,0,1)[12] - winner of the {p,q:0-2, P,Q:0-1} grid at
# d=1, D=0 (05_sarima.R). fitdf = p+q+P+Q = 4.
fits[["SARIMA"]] <- Arima(train, order = c(0, 1, 2),
  seasonal = list(order = c(1, 0, 1), period = 12))

# ETS(A,N,A) - winner of the 9-candidate additive/multiplicative x
# damped/undamped search (06_ets.R). No trend, so damped is moot (NULL).
fits[["ETS"]] <- ets(train, model = "ANA", damped = NULL)

# BATS(no box-cox, damped trend) - winner of the 4-config
# {box-cox, damped} x {on, off} grid (07_bats.R).
fits[["BATS"]] <- bats(train, use.box.cox = FALSE, use.trend = TRUE,
                       use.damped.trend = TRUE)

# ---- Build the comparison table ---------------------------------------
rows <- data.frame()
for (nm in names(fits)) {
  fit <- fits[[nm]]
  fc  <- if (nm == "SNAIVE") fit else if (nm == "ARIMA") fc_arima_precomputed else forecast(fit, h = H)
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
  if (nm == "SNAIVE") fits[[nm]]
  else if (nm == "ARIMA") fc_arima_precomputed
  else forecast(fits[[nm]], h = H))
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

# ---- 2x2 panel: each family member's own forecast vs actual -----------
# Same idea as the combined overlay plot above, but one panel per model
# instead of five lines on one axis - easier to read each model's own
# fit/miss pattern without the others crowding it. Same grid.arrange()
# (preview) + arrangeGrob() (save) pattern as 03_eda.R's 2x2 EDA figure.
zoom_start <- c(2022, 1)   # matches the combined overlay plot above

panel <- function(nm) {
  autoplot(window(y, start = zoom_start)) +
    autolayer(fc_all[[nm]]$mean, series = "Forecast", size = 0.8) +
    autolayer(test, series = "Actual", size = 0.8) +
    scale_colour_manual(values = c(Forecast = "steelblue", Actual = "red")) +
    ggtitle(nm) +
    ylab("Monthly ridership") + scale_y_millions() +
    theme(legend.position = "none", plot.title = element_text(size = 11))
}

p_arima_panel  <- panel("ARIMA")
p_sarima_panel <- panel("SARIMA")
p_ets_panel    <- panel("ETS")
p_bats_panel   <- panel("BATS")

grid.arrange(p_arima_panel, p_sarima_panel, p_ets_panel, p_bats_panel, ncol = 2,
             top = "Forecast vs Actual by Model (12-month holdout)")

p_grid_2x2 <- arrangeGrob(p_arima_panel, p_sarima_panel, p_ets_panel, p_bats_panel,
                           ncol = 2, top = "Forecast vs Actual by Model (12-month holdout)")
ggsave(fig("comparison", "model_grid_2x2.png"), p_grid_2x2, width = 11, height = 7, dpi = 150)
