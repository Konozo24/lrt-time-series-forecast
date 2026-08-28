# 00_setup.R - packages. Run once per Posit Cloud session.
# Free tier: 1GB RAM. Do NOT install prophet (Stan compile fails/times out).

pkgs <- c("forecast", "tseries", "dplyr", "ggplot2", "lubridate", "gridExtra")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

library(forecast)   # tbats(), accuracy(), checkresiduals()
library(tseries)    # adf.test(), kpss.test()
library(dplyr)
library(ggplot2)
library(lubridate)
library(gridExtra)  # grid.arrange() - combine multiple ggplot panels into one figure

# Topic: LRT Ampang line monthly ridership, SDG 11 (Sustainable Cities and
# Communities), Target 11.2. Full period 2019-01 to 2026-06 retained per
# tutor's instruction (do NOT truncate out the MCO period) - see
# 02_mco_resolution.R for how the 2020-2021 disruption is handled instead
# of removed.
