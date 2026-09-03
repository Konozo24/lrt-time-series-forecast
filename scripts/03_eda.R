# 03_eda.R - EDA + stationarity/seasonality checks on the full-period
# (55-month) Rapid KL bus series. No reconstruction step precedes this -
# unlike an earlier version of this project, this series has no COVID-era
# gap to resolve (see 01_data_prep.R).

source("scripts/00_setup.R")
bus_ts <- readRDS("data/bus_rkl_monthly.rds")

# --- Visual EDA: combined into ONE figure (2x2 grid) ---
# Panel 1: raw series, Panel 2: seasonal plot, Panel 3: ACF, Panel 4: PACF
# (all on the differenced series for ACF/PACF, since that's what informs
# ARIMA/SARIMA order selection.

p1 <- autoplot(bus_ts) +
  ggtitle("Monthly Rapid KL Bus Ridership (2022-2026)") +
  theme(plot.title = element_text(size = 10)) + scale_y_millions()

p2 <- ggseasonplot(bus_ts) +
  ggtitle("Seasonal Plot") +
  theme(plot.title = element_text(size = 10), legend.position = "none") +
  scale_y_millions()

p3 <- ggAcf(diff(bus_ts), lag.max = 24) +
  ggtitle("ACF (first difference)") +
  theme(plot.title = element_text(size = 10))

p4 <- ggPacf(diff(bus_ts), lag.max = 24) +
  ggtitle("PACF (first difference)") +
  theme(plot.title = element_text(size = 10))

grid.arrange(p1, p2, p3, p4, ncol = 2,
             top = "Rapid KL Bus - EDA Summary")

# Save the combined figure directly (useful for pasting into the report)
combined <- arrangeGrob(p1, p2, p3, p4, ncol = 2, top = "Rapid KL Bus - EDA Summary")
ggsave(fig("eda", "eda_summary.png"), combined, width = 10, height = 7, dpi = 150)

# STL decomposition
decomp <- stl(bus_ts, s.window = "periodic", robust = TRUE)
print(autoplot(decomp) + ggtitle("STL Decomposition") + scale_y_millions())
ggsave(fig("eda", "stl_decomposition.png"), width = 8, height = 6, dpi = 150)

# Seasonal plot, also saved standalone (p2 above is one quadrant of the
# 2x2 grid - too small to read on its own in the report; this is the
# same plot at full size)
ggsave(fig("eda", "seasonal_plot.png"), p2, width = 8, height = 6, dpi = 150)

# Lag plot (Hyndman-style): scatterplots of the series against itself at
# lags 1-12. A strong positive linear pattern specifically at lag 12
# (relative to the others) is the visual counterpart to the ACF/PACF
# lag-12 spike above - confirms the annual seasonality directly as a
# shape, not just a bar height. Kept as its own figure (also already a
# multi-panel grid on its own, same reasoning as STL above).
print(gglagplot(bus_ts, lags = 12) + ggtitle("Lag Plot (lags 1-12)"))
ggsave(fig("eda", "lag_plot.png"), width = 8, height = 8, dpi = 150)

# Seasonal strength (Hyndman's measure)
seasonal_strength <- max(0, 1 - var(decomp$time.series[, "remainder"]) /
                            var(decomp$time.series[, "seasonal"] + decomp$time.series[, "remainder"]))
cat("Seasonal strength:", round(seasonal_strength, 3), "\n")

# --- Stationarity tests ---
cat("\n--- Level series ---\n")
print(adf.test(bus_ts))   # want p < 0.05 for stationary
print(kpss.test(bus_ts))  # want p > 0.05 for stationary

cat("\n--- First difference ---\n")
print(adf.test(diff(bus_ts)))
print(kpss.test(diff(bus_ts)))
