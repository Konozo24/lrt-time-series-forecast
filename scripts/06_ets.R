pkgs <- c("forecast", "tseries", "dplyr", "lubridate", "ggplot2", "zoo", "Kendall", "arrow")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(arrow)      
library(forecast)   
library(tseries)    
library(dplyr)
library(lubridate)
library(ggplot2)

dir.create("output/plots/ets", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables",    showWarnings = FALSE, recursive = TRUE)
dir.create("output/models",    showWarnings = FALSE, recursive = TRUE)

fig <- function(category, name) file.path("output/plots", category, name)
tbl <- function(name) file.path("output/tables", name)
mdl <- function(name) file.path("output/models", name)

scale_y_millions <- function() {
  scale_y_continuous(labels = function(x) paste0(round(x / 1e6, 1), "M"))
}

# ---------------------------------------------------------------------
raw <- as.data.frame(read_parquet("https://storage.data.gov.my/transportation/ridership_headline.parquet"))
raw$date <- as.Date(raw$date)

monthly <- raw %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(rail_lrt_ampang = sum(rail_lrt_ampang, na.rm = TRUE),
            n_days = n()) %>%
  ungroup() %>%
  filter(n_days >= 28)   # drop partial calendar months

cat("Monthly series:", nrow(monthly), "months,",
    format(min(monthly$month)), "to", format(max(monthly$month)), "\n")

start_year  <- year(min(monthly$month))
start_month <- month(min(monthly$month))
ampang_ts   <- ts(monthly$rail_lrt_ampang, start = c(start_year, start_month), frequency = 12)

# ---------------------------------------------------------------------
mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")
mco_mask  <- monthly$month >= mco_start & monthly$month <= mco_end
cat("MCO window:", sum(mco_mask), "months (", format(mco_start), "to", format(mco_end), ")\n")

temp_filled <- monthly$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)
stl_temp     <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est <- as.numeric(stl_temp$time.series[, "seasonal"])

i_pre  <- which(monthly$month == mco_start - months(1))
i_post <- which(monthly$month == mco_end + months(1))
last_pre   <- monthly$rail_lrt_ampang[i_pre]  - seasonal_est[i_pre]     # seasonally adjusted
first_post <- monthly$rail_lrt_ampang[i_post] - seasonal_est[i_post]
n_gap <- sum(mco_mask)
bridge_trend <- seq(last_pre, first_post, length.out = n_gap + 2)[2:(n_gap + 1)]

resolved <- monthly$rail_lrt_ampang
resolved[mco_mask] <- bridge_trend + seasonal_est[mco_mask]

cat("Bridged", n_gap, "MCO-window months (trend interpolation + STL seasonal),",
    "0 months dropped.\n")

ampang_resolved <- ts(resolved, start = start(ampang_ts), frequency = 12)

p_mco <- autoplot(cbind(Original = ampang_ts, Resolved = ampang_resolved)) +
  ggtitle("LRT Ampang: Original vs. MCO-Resolved Series") +
  ylab("Monthly ridership") + scale_y_millions()
print(p_mco)
ggsave(fig("ets", "mco_resolved.png"), p_mco, width = 8, height = 5, dpi = 150)

y <- ampang_resolved   # series used for everything from here on

# ---------------------------------------------------------------------
p_stl <- autoplot(stl(y, s.window = "periodic", robust = TRUE)) +
  ggtitle("STL decomposition")
print(p_stl)
ggsave(fig("ets", "stl_decomposition.png"), p_stl, width = 8, height = 5, dpi = 150)

cat("\n== ADF (want p < 0.05 for stationary) ==\n")
print(adf.test(y))

cat("\n== KPSS (want p > 0.05 for stationary) ==\n")
print(kpss.test(y))

cat("\n== Ljung-Box on RAW series (want p < 0.05 -> not white noise) ==\n")
print(Box.test(y, lag = 12, type = "Ljung-Box"))
print(Box.test(y, lag = 24, type = "Ljung-Box"))

if (requireNamespace("Kendall", quietly = TRUE)) {
  cat("\n== Mann-Kendall trend test (H0: no monotonic trend) ==\n")
  print(Kendall::MannKendall(as.numeric(y)))
  cat("tau near 0 (|tau| < 0.1) = negligible trend magnitude even if p is\n",
      "small; 0.1-0.3 = weak.\n")
}

cat("\n== Min value (near-zero check for MAPE stability) ==\n")
print(min(y, na.rm = TRUE))

# ---------------------------------------------------------------------
h     <- 12
H     <- h   # alias - the ETS candidate loop below uses H
train <- head(y, length(y) - h)
test  <- tail(y, h)

LAG_MAX <- 16

acf_out_of_bounds <- function(resid, lag.max = LAG_MAX) {
  r  <- na.omit(resid)
  n  <- length(r)
  ci <- 1.96 / sqrt(n)
  a  <- acf(r, plot = FALSE, lag.max = lag.max)$acf[-1]
  sum(abs(a) > ci)
}

model_fitdf <- function(fit) {
  if (!is.null(fit$arma) && length(fit$arma) >= 4) sum(fit$arma[1:4]) else 0
}
lb_test <- function(fit, lag) {
  Box.test(residuals(fit), lag = lag, fitdf = model_fitdf(fit), type = "Ljung-Box")
}

gap_check <- function(acc_matrix) {
  mape_tr <- acc_matrix["Training set", "MAPE"]
  mape_te <- acc_matrix["Test set", "MAPE"]
  mase_tr <- acc_matrix["Training set", "MASE"]
  mase_te <- acc_matrix["Test set", "MASE"]
  rmse_tr <- acc_matrix["Training set", "RMSE"]
  rmse_te <- acc_matrix["Test set", "RMSE"]
  mape_gap   <- abs(mape_te - mape_tr) / mape_te
  mase_gap   <- abs(mase_te - mase_tr) / mase_te
  rmse_ratio <- rmse_te / rmse_tr
  list(mape_gap = mape_gap, mase_gap = mase_gap, rmse_ratio = rmse_ratio,
       rmse_train = rmse_tr, mase_train = mase_tr,
       direction = if (mase_te > mase_tr) "test worse" else "train worse",
       within_10pct = (mase_gap <= 0.10), within_1_3x = (rmse_ratio <= 1.3))
}

model_summary <- function(name, fit, fc, test) {
  a  <- accuracy(fc, test)
  g  <- gap_check(a)
  r  <- residuals(fit)
  lb12  <- lb_test(fit, 12)
  lbmax <- lb_test(fit, LAG_MAX)
  data.frame(
    model = name,
    MAPE_test = round(a["Test set", "MAPE"], 3),
    RMSE_test = round(a["Test set", "RMSE"], 0),
    MAE_test  = round(a["Test set", "MAE"], 0),
    MASE_test = round(a["Test set", "MASE"], 3),
    RMSE_train = round(g$rmse_train, 0),
    MASE_train = round(g$mase_train, 3),
    fitdf = model_fitdf(fit),
    lb_pvalue_12 = round(lb12$p.value, 4),
    lb_pvalue_16 = round(lbmax$p.value, 4),
    lb_pass = (lb12$p.value > 0.05) & (lbmax$p.value > 0.05),
    n_lags_out_12 = acf_out_of_bounds(r, lag.max = 12),
    n_lags_out_16 = acf_out_of_bounds(r, lag.max = LAG_MAX),
    mase_gap_pct = round(g$mase_gap * 100, 1),
    within_10pct = g$within_10pct,
    rmse_ratio = round(g$rmse_ratio, 3),
    within_1_3x = g$within_1_3x,
    direction = g$direction,
    mape_gap_pct = round(g$mape_gap * 100, 1),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------
candidates <- list(
  list(label = "ETS(A,N,N)",     model = "ANN", damped = NULL,  lambda = NULL),
  list(label = "ETS(A,A,N)",     model = "AAN", damped = FALSE, lambda = NULL),
  list(label = "ETS(A,Ad,N)",    model = "AAN", damped = TRUE,  lambda = NULL),
  list(label = "ETS(A,N,A)",     model = "ANA", damped = NULL,  lambda = NULL),
  list(label = "ETS(A,A,A)",     model = "AAA", damped = FALSE, lambda = NULL),
  list(label = "ETS(A,Ad,A)",    model = "AAA", damped = TRUE,  lambda = NULL),
  list(label = "ETS(M,N,M)",     model = "MNM", damped = NULL,  lambda = NULL),
  list(label = "ETS(M,A,M)",     model = "MAM", damped = FALSE, lambda = NULL),
  list(label = "ETS(M,Ad,M)",    model = "MAM", damped = TRUE,  lambda = NULL),
  # Box-Cox variants to bring residual spikes inside blue dashed lines
  list(label = "ETS(A,A,A)+BC",  model = "AAA", damped = FALSE, lambda = "auto"),
  list(label = "ETS(M,A,M)+BC",  model = "MAM", damped = FALSE, lambda = "auto")
)

results <- data.frame()
fits <- list()
forecasts <- list()

for (cand in candidates) {
  fit <- tryCatch(
    ets(train, model = cand$model, damped = cand$damped, lambda = cand$lambda),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    fc  <- forecast(fit, h = H)
    acc <- accuracy(fc, test)
    
    train_rmse <- acc["Training set", "RMSE"]
    test_rmse  <- acc["Test set", "RMSE"]
    test_mae   <- acc["Test set", "MAE"]
    test_mape  <- acc["Test set", "MAPE"]
    test_mase  <- acc["Test set", "MASE"]
    
    rmse_ratio <- test_rmse / train_rmse
    
    lb_res    <- Box.test(residuals(fit), lag = LAG_MAX, type = "Ljung-Box")
    lb_pvalue <- lb_res$p.value
    
    res <- na.omit(residuals(fit))
    res_acf <- acf(res, lag.max = LAG_MAX, plot = FALSE)$acf[-1]
    acf_threshold <- 1.96 / sqrt(length(res))
    acf_out_count <- sum(abs(res_acf) > acf_threshold)
    
    meets_criteria <- (rmse_ratio < 1.3) &&
      (test_mape < 10) &&
      (test_mase < 1) &&
      (!is.na(lb_pvalue) && lb_pvalue > 0.05) &&
      (acf_out_count <= 2)
    
    fits[[cand$label]]      <- fit
    forecasts[[cand$label]] <- fc
    
    results <- rbind(results, data.frame(
      spec          = cand$label,
      AICc          = fit$aicc,
      BIC           = fit$bic,
      RMSE_Ratio    = round(rmse_ratio, 3),
      RMSE          = round(test_rmse, 2),
      MAE           = round(test_mae, 2),
      MAPE          = round(test_mape, 2),
      MASE          = round(test_mase, 3),
      LB_pvalue     = round(lb_pvalue, 4),
      ACF_Out_Lags  = acf_out_count,
      Valid         = meets_criteria
    ))
  }
}

cat("\n--- All ETS Candidate Models Evaluated ---\n")
print(results, row.names = FALSE)

if (nrow(results) == 0) {
  stop("Error: No models were successfully fitted. Check your training data or candidates list.")
}

valid_results <- results[results$Valid == TRUE, ]

if (nrow(valid_results) > 0) {
  valid_results <- valid_results[order(valid_results$AICc), ]
  best_label <- valid_results$spec[1]
  cat("\nSelected Best Valid Model:", best_label, " (AICc =", round(valid_results$AICc[1], 2), ")\n")
} else {
  warning("No ETS model satisfied ALL constraints. Falling back to model with lowest AICc.")
  results_sorted <- results[order(results$AICc), ]
  best_label <- results_sorted$spec[1]
  cat("\nSelected Fallback Model:", best_label, "\n")
}

fit_ets <- fits[[best_label]]
print(summary(fit_ets))

# Zoom parameters (last 3 years)
end_year   <- end(test)[1]
start_zoom <- c(end_year - 3, end(test)[2])
train_zoom <- window(train, start = start_zoom)

fc_ets <- forecast(fit_ets, h = H)

p_fc <- autoplot(train_zoom, series = "Training Data") +
  autolayer(test, series = "Actual (Test Data)") +
  autolayer(fc_ets$mean, series = "ETS Forecast") +
  ggtitle(paste(best_label, "Forecast vs Actual (Zoomed: Last 3 Years)"))

print(p_fc)
ggsave(fig("ets", "ets_forecast.png"), p_fc, width = 8, height = 5, dpi = 150)

cat("\n--- Holdout accuracy (12-month full seasonal cycle) ---\n")
print(accuracy(fc_ets, test))

cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_ets)
dev.copy(png, filename = fig("ets", "ets_residuals.png"), width = 900, height = 700, res = 130)
dev.off()

print(Box.test(residuals(fit_ets), lag = LAG_MAX, type = "Ljung-Box"))
cat("ACF-out-of-bounds:", acf_out_of_bounds(residuals(fit_ets)), "/", LAG_MAX, "\n")

g <- gap_check(accuracy(fc_ets, test))
cat("Train/test MAPE gap:", round(g$mape_gap * 100, 1),
    "% | within 10%:", g$within_10pct, "\n")

summ <- model_summary(best_label, fit_ets, fc_ets, test)
print(summ)

write.csv(results, tbl("ets_grid.csv"), row.names = FALSE)
write.csv(summ, tbl("summary_ets.csv"), row.names = FALSE)
saveRDS(fit_ets, mdl("fit_ets.rds"))
saveRDS(fc_ets,  mdl("fc_ets.rds"))
