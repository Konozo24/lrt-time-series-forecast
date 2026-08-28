# 02_mco_resolution.R - replace the COVID-19 MCO months with an estimate of
# what ridership would have been, keeping all months in the series.
#
# For the disrupted window (Mar 2020 - Dec 2021):
#   trend    = linear interpolation between Feb 2020 and Jan 2022
#   seasonal = STL seasonal component of the series
#   resolved = trend + seasonal

source("scripts/00_setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full.rds")
monthly_df <- readRDS("data/ampang_monthly_full_df.rds")

mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")
mco_mask  <- monthly_df$month >= mco_start & monthly_df$month <= mco_end
cat("MCO window:", sum(mco_mask), "months (", format(mco_start), "to", format(mco_end), ")\n")

# Linear-fill the window so STL has a gap-free series; only the seasonal
# component of this fit is used below
temp_filled <- monthly_df$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)
stl_temp <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est <- as.numeric(stl_temp$time.series[, "seasonal"])

# Interpolate the trend between the last pre-MCO and first post-MCO values
last_pre  <- monthly_df$rail_lrt_ampang[monthly_df$month == mco_start - months(1)]
first_post <- monthly_df$rail_lrt_ampang[monthly_df$month == mco_end + months(1)]
n_gap <- sum(mco_mask)
bridge_trend <- seq(last_pre, first_post, length.out = n_gap + 2)[2:(n_gap + 1)]

cat("Bridging trend from", format(mco_start - months(1)), "(", last_pre, ") to",
    format(mco_end + months(1)), "(", first_post, ")\n")

resolved <- monthly_df$rail_lrt_ampang
resolved[mco_mask] <- bridge_trend + seasonal_est[mco_mask]

comparison <- data.frame(
  month = monthly_df$month[mco_mask],
  original = monthly_df$rail_lrt_ampang[mco_mask],
  resolved = round(resolved[mco_mask])
)
print(comparison)

ampang_ts_resolved <- ts(resolved, start = start(ampang_ts), frequency = 12)
saveRDS(ampang_ts_resolved, "data/ampang_monthly_full_resolved.rds")

cat("\nSaved data/ampang_monthly_full_resolved.rds - length:", length(ampang_ts_resolved),
    "(all", length(ampang_ts), "months retained, MCO window resolved not removed)\n")

print(
  autoplot(cbind(Original = ampang_ts, Resolved = ampang_ts_resolved)) +
    ggtitle("LRT Ampang: Original vs. MCO-Resolved Series") +
    ylab("Monthly ridership")
)
