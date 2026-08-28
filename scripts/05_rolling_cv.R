# 05_rolling_cv.R - rolling-origin CV for TBATS, since a single 6-month
# holdout can land on an atypical window. Expanding training window,
# 6-month forecast horizon, 7 folds (mirrors the group's shared CV
# methodology).

source("scripts/00_setup.R")
ampang_ts <- readRDS("data/ampang_monthly_full_resolved.rds")

h <- 6
n <- length(ampang_ts)
train_sizes <- seq(n - h - 18, n - h, by = 3)  # 7 folds ending at n-h

mapes <- c()
for (ts_size in train_sizes) {
  tr <- head(ampang_ts, ts_size)
  te <- window(ampang_ts, start = time(ampang_ts)[ts_size + 1],
               end = time(ampang_ts)[ts_size + h])
  fit_cv <- tbats(tr, use.box.cox = TRUE, use.trend = TRUE, use.damped.trend = TRUE)
  fc_cv <- forecast(fit_cv, h = h)
  mape <- mean(abs((te - fc_cv$mean) / te)) * 100
  mapes <- c(mapes, mape)
  cat("train size:", ts_size, " MAPE:", round(mape, 2), "%\n")
}

cat("\nRolling-CV mean MAPE:", round(mean(mapes), 2), "%  std:", round(sd(mapes), 2), "\n")

write.csv(data.frame(train_size = train_sizes, mape = mapes),
          "output/tbats_rolling_cv.csv", row.names = FALSE)
