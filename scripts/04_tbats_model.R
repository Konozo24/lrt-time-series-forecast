# 04_tbats_model.R - Ming's model family: ETS (TBATS variant).
#
# Model pick: TBATS(box-cox, trend, damped trend, seasonal period 12).
# Selected based on two EDA-derived hypotheses (see project report /
# 03_eda.R output for supporting numbers):
#   1. A multiplicative-vs-additive seasonal ETS comparison found
#      multiplicative seasonality fits better - the seasonal swing scales
#      proportionally with ridership level. TBATS's Box-Cox transform is
#      designed to address exactly this kind of variance pattern.
#   2. TBATS represents seasonality with a small number of Fourier
#      (trigonometric) terms rather than 12 separate monthly seasonal-
#      state parameters, which is more parameter-efficient - relevant
#      given the dataset's sample size.
#
# Fitted on the FULL 91-month, MCO-resolved series (see 02_mco_resolution.R)
# per tutor's instruction that all periods must be included.

source("scripts/00_setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full_resolved.rds")

h <- 6
train <- head(ampang_ts, length(ampang_ts) - h)
test  <- tail(ampang_ts, h)

fit_tbats <- tbats(train, use.box.cox = TRUE, use.trend = TRUE, use.damped.trend = TRUE)
print(fit_tbats)

fc_tbats <- forecast(fit_tbats, h = h)
print(autoplot(fc_tbats) + autolayer(test, series = "Actual"))

cat("\n--- Holdout accuracy (last 6 months) ---\n")
print(accuracy(fc_tbats, test))

cat("\n--- Residual diagnostics ---\n")
checkresiduals(fit_tbats)
print(Box.test(residuals(fit_tbats), lag = 6, type = "Ljung-Box"))
print(Box.test(residuals(fit_tbats), lag = 12, type = "Ljung-Box"))

saveRDS(fit_tbats, "output/fit_tbats.rds")
saveRDS(fc_tbats, "output/fc_tbats.rds")

# Expected reference numbers (from Python prototype, confirm your R
# numbers land close to these): holdout MAPE ~4.97%, Ljung-Box p~0.36
# (lag 6) / ~0.38 (lag 12) - both clear passes on the full resolved series,
# notably better than fitting only on a truncated 2022+ subset (MAPE
# ~7.24%, borderline Ljung-Box at lag 12) - more retained data lets TBATS
# estimate trend/seasonal structure more reliably.
