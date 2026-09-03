# 01_data_prep.R - load + aggregate Rapid KL bus daily ridership to monthly.
# Source: data.gov.my Daily Public Transport Ridership dataset
# (https://data.gov.my/data-catalogue/ridership_headline?visual=bus_rkl),
# Prasarana Malaysia + Ministry of Transport, CC BY 4.0.
#
# Column: bus_rkl (Rapid KL bus). This column has no data before 2022-01,
# so the series is filtered to 2022-01 onward explicitly - this is a
# genuinely complete series from its first observation, not a truncation
# of a longer one. There is no COVID-era gap to resolve: unlike the rail
# series used in an earlier version of this project, bus_rkl simply does
# not report anything for the disrupted period, so no reconstruction step
# is needed (there is no 02_mco_resolution.R in this version).

source("scripts/00_setup.R")

raw <- read.csv("data/ridership_headline.csv", stringsAsFactors = FALSE)
raw$date <- as.Date(raw$date)
raw <- raw[raw$date >= as.Date("2022-01-01"), ]

# Aggregate daily -> monthly totals (sum of trips per calendar month)
monthly <- raw %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    bus_rkl = sum(bus_rkl, na.rm = TRUE),
    n_days = n()
  ) %>%
  ungroup()

# Drop partial calendar month (e.g. a mid-month data pull)
monthly <- monthly %>% filter(n_days >= 28)

cat("Monthly series:", nrow(monthly), "months,",
    format(min(monthly$month)), "to", format(max(monthly$month)), "\n")

# Convert to ts object, frequency = 12 (monthly seasonality)
start_year  <- year(min(monthly$month))
start_month <- month(min(monthly$month))
bus_ts <- ts(monthly$bus_rkl, start = c(start_year, start_month), frequency = 12)

saveRDS(bus_ts,  "data/bus_rkl_monthly.rds")
saveRDS(monthly, "data/bus_rkl_monthly_df.rds")

cat("Saved data/bus_rkl_monthly.rds - length:", length(bus_ts), "\n")
