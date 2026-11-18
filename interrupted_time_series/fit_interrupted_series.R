#!/usr/bin/env Rscript
# Interrupted time series analysis of annual resistance prevalence.
# The intervention is announcement 2292, which took effect at the end of
# 2016, so 2017 is the first post-intervention year.

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
out_table <- args[2]
out_figure <- args[3]
intervention_year <- if (length(args) >= 4) as.integer(args[4]) else 2017

data <- read.delim(input, sep = "\t", check.names = FALSE)
data <- data[order(data$year), ]
data$time <- seq_len(nrow(data)) - 1
data$intervention <- ifelse(data$year >= intervention_year, 1, 0)
data$time_after <- ifelse(data$intervention == 1,
                          data$year - intervention_year + 1, 0)

determinants <- c("qrdr", "gyrA", "parC", "qnrS1", "oqxAB")
results <- list()

for (det in determinants) {
  model_data <- data.frame(
    prevalence = data[[paste0(det, "_prevalence")]],
    time = data$time,
    intervention = data$intervention,
    time_after = data$time_after,
    total = data$total)
  fit <- lm(prevalence ~ time + intervention + time_after,
            data = model_data, weights = model_data$total)
  coef_table <- coeftest(fit, vcov = NeweyWest(fit, lag = 1,
                                              prewhite = FALSE))
  coefficients <- coef(fit)
  results[[det]] <- data.frame(
    determinant = det,
    pre_slope = coefficients["time"],
    level_change = coefficients["intervention"],
    slope_change = coefficients["time_after"],
    post_slope = coefficients["time"] + coefficients["time_after"],
    level_change_p = coef_table["intervention", "Pr(>|t|)"],
    slope_change_p = coef_table["time_after", "Pr(>|t|)"],
    row.names = NULL)
  data[[paste0(det, "_fitted")]] <- fitted(fit)
}

summary_table <- do.call(rbind, results)
write.table(summary_table, out_table, sep = "\t", quote = FALSE,
            row.names = FALSE)
print(summary_table, digits = 3)

png(out_figure, width = 1800, height = 1200, res = 200)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (det in determinants) {
  observed <- data[[paste0(det, "_prevalence")]]
  fitted <- data[[paste0(det, "_fitted")]]
  plot(data$year, observed, pch = 19, xlab = "Year",
       ylab = "Prevalence", main = det, ylim = c(0, max(observed, 1)))
  lines(data$year, fitted, col = "#4575B4", lwd = 2)
  abline(v = intervention_year - 0.5, lty = 2, col = "#e54e21")
}
plot.new()
dev.off()
