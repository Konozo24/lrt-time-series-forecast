# 03_eda.R - EDA + stationarity/seasonality checks on the MCO-resolved,
# full-period (91-month) series.

source("scripts/00_setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full_resolved.rds")

# --- Visual EDA ---
autoplot(ampang_ts) + ggtitle("LRT Ampang Monthly Ridership (2019-2026, MCO-resolved)")
ggseasonplot(ampang_ts)

decomp <- stl(ampang_ts, s.window = "periodic", robust = TRUE)
autoplot(decomp)

# Seasonal strength (Hyndman's measure)
seasonal_strength <- max(0, 1 - var(decomp$time.series[, "remainder"]) /
                            var(decomp$time.series[, "seasonal"] + decomp$time.series[, "remainder"]))
cat("Seasonal strength:", round(seasonal_strength, 3), "\n")

# --- Stationarity tests ---
cat("\n--- Level series ---\n")
print(adf.test(ampang_ts))   # want p < 0.05 for stationary
print(kpss.test(ampang_ts))  # want p > 0.05 for stationary

cat("\n--- First difference ---\n")
print(adf.test(diff(ampang_ts)))
print(kpss.test(diff(ampang_ts)))

# --- ACF / PACF ---
ggAcf(diff(ampang_ts), lag.max = 24)
ggPacf(diff(ampang_ts), lag.max = 24)
