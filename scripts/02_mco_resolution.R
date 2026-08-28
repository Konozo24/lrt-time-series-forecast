# 02_mco_resolution.R - resolve the COVID-19 MCO structural break WITHOUT
# removing any months from the series (tutor's requirement: all periods,
# including MCO, must be included and the resolution explained).
#
# Approach tried FIRST and REJECTED: automatic outlier detection via STL
# residual thresholding (flag any month with |remainder| > 3*MAD). This
# produced false positives - it incorrectly flagged Jan-Feb 2020 (genuine
# normal pre-lockdown ridership, ~5.1-5.4M) and June 2022 (genuine recovery
# month) as anomalies, because a single sharp discontinuity distorts STL's
# local trend estimate right at its edges, not just during the disruption
# itself. Automatic residual-threshold detection is therefore NOT robust
# to one large structural break - it's built for scattered, independent
# outliers, not a multi-month regime shift.
#
# Approach used INSTEAD: a known-intervention adjustment based on Malaysia's
# actual MCO timeline (domain knowledge), not statistical detection:
#   1. Explicitly define the disrupted window: March 2020 - December 2021
#      (MCO 1.0 through the National Recovery Plan / endemic transition),
#      22 months.
#   2. Bridge the TREND across this window with linear interpolation
#      between the last stable pre-MCO value (Feb 2020) and the first
#      stable post-recovery value (Jan 2022) - this avoids using the
#      collapsed values themselves to estimate what the trend "should"
#      have looked like.
#   3. Add back a SEASONAL component estimated from an STL fit on the
#      non-disrupted months only (with the disrupted window temporarily
#      linearly filled just so STL has no gaps to fit through).
#   4. Replace only the 22 disrupted months with (bridged trend + seasonal
#      estimate). All 91 months remain in the series - none are dropped -
#      but the collapse-and-recovery no longer distorts trend/seasonal
#      estimation for the models fitted on this series.
#
# This is a defensible, explainable resolution: it keeps the full period
# the tutor requires, and the choice of window is grounded in known
# real-world dates rather than an automatic statistical rule that we
# verified (above) is not reliable for this kind of single, sharp shock.

source("scripts/00_setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full.rds")
monthly_df <- readRDS("data/ampang_monthly_full_df.rds")

mco_start <- as.Date("2020-03-01")
mco_end   <- as.Date("2021-12-01")
mco_mask  <- monthly_df$month >= mco_start & monthly_df$month <= mco_end
cat("MCO window:", sum(mco_mask), "months (", format(mco_start), "to", format(mco_end), ")\n")

# Temporarily linear-fill the disrupted window just so STL has a complete
# series to extract a seasonal estimate from (the filled values here are
# NOT the final resolution - only the seasonal component from this fit is
# reused below)
temp_filled <- monthly_df$rail_lrt_ampang
temp_filled[mco_mask] <- NA
temp_filled_ts <- ts(zoo::na.approx(temp_filled), start = start(ampang_ts), frequency = 12)
stl_temp <- stl(temp_filled_ts, s.window = "periodic", robust = TRUE)
seasonal_est <- as.numeric(stl_temp$time.series[, "seasonal"])

# Bridge the trend linearly between the last pre-MCO point and first
# post-MCO point (does not use any collapsed values as anchors)
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
