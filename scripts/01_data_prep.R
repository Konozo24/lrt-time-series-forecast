# 01_data_prep.R - load + aggregate LRT Ampang daily ridership to monthly.
# Source: data.gov.my Daily Public Transport Ridership dataset
# (https://data.gov.my/data-catalogue/ridership_headline), Prasarana
# Malaysia + Ministry of Transport, CC BY 4.0.
#
# IMPORTANT: full period retained (2019-01 to 2026-06), including the
# COVID-19 MCO period, per tutor's instruction. The 2020-2021 disruption
# is resolved via a known-intervention adjustment in 02_mco_resolution.R,
# not by truncating the series.

source("scripts/00_setup.R")

raw <- read.csv("data/ridership_headline.csv", stringsAsFactors = FALSE)
raw$date <- as.Date(raw$date)

# Aggregate daily -> monthly totals (sum of trips per calendar month)
monthly <- raw %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    rail_lrt_ampang = sum(rail_lrt_ampang, na.rm = TRUE),
    n_days = n()
  ) %>%
  ungroup()

# Drop the last month if it's incomplete (partial calendar month at the
# time of data extraction)
monthly <- monthly %>% filter(n_days >= 28)

cat("Monthly series:", nrow(monthly), "months,",
    format(min(monthly$month)), "to", format(max(monthly$month)), "\n")

# Convert to ts object, frequency = 12 (monthly seasonality)
start_year  <- year(min(monthly$month))
start_month <- month(min(monthly$month))
ampang_ts <- ts(monthly$rail_lrt_ampang, start = c(start_year, start_month), frequency = 12)

saveRDS(ampang_ts, "data/ampang_monthly_full.rds")
saveRDS(monthly, "data/ampang_monthly_full_df.rds")

cat("Saved data/ampang_monthly_full.rds - length:", length(ampang_ts), "\n")
