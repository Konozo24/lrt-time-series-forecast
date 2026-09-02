# 03_eda.R - EDA + stationarity/seasonality checks on the MCO-resolved,
# full-period (91-month) series.

source("scripts/00_setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full_resolved.rds")

# --- Visual EDA: combined into ONE figure (2x2 grid) ---
# Panel 1: raw series, Panel 2: seasonal plot, Panel 3: ACF, Panel 4: PACF
# (all on the differenced series for ACF/PACF, since that's what informs
# ARIMA/SARIMA order selection.

p1 <- autoplot(ampang_ts) +
  ggtitle("Monthly Ridership (2019-2026, MCO-resolved)") +
  theme(plot.title = element_text(size = 10))

p2 <- ggseasonplot(ampang_ts) +
  ggtitle("Seasonal Plot") +
  theme(plot.title = element_text(size = 10), legend.position = "none")

p3 <- ggAcf(diff(ampang_ts), lag.max = 24) +
  ggtitle("ACF (first difference)") +
  theme(plot.title = element_text(size = 10))

p4 <- ggPacf(diff(ampang_ts), lag.max = 24) +
  ggtitle("PACF (first difference)") +
  theme(plot.title = element_text(size = 10))

grid.arrange(p1, p2, p3, p4, ncol = 2,
             top = "LRT Ampang - EDA Summary")

# Save the combined figure directly (useful for pasting into the report)
combined <- arrangeGrob(p1, p2, p3, p4, ncol = 2, top = "LRT Ampang - EDA Summary")
ggsave(fig("eda", "eda_summary.png"), combined, width = 10, height = 7, dpi = 150)

# STL decomposition
decomp <- stl(ampang_ts, s.window = "periodic", robust = TRUE)
print(autoplot(decomp) + ggtitle("STL Decomposition"))
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
print(gglagplot(ampang_ts, lags = 12) + ggtitle("Lag Plot (lags 1-12)"))
ggsave(fig("eda", "lag_plot.png"), width = 8, height = 8, dpi = 150)

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
