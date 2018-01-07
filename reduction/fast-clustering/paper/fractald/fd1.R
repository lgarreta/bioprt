#!/usr/bin/Rscript

## Not run:
library(RandomFields)
n <- 10000

# generate a time series
rf <- GaussRF(x = c(0, 1, 1/n), model = "stable",
              grid = TRUE, gridtriple = TRUE,
              param = c(mean=0, variance=1, nugget=0, scale=100, kappa=1))

# Plots for two sliding windows of each of the four methods below.
# Argument nlags is common to all methods;
# the 'variation' method has in addition argument p.index
par(mfrow=c(2,4)) # one row per window
fd <- fd.estimate(rf,
                  methods = list(list(name="variation", p.index=0.5),
                  "variogram", "hallwood", "boxcount"),
                  window.size = 5000, step.size = 5000, plot.loglog = TRUE, nlags = 10)

# 2d random fields
n <- 200
rf2d <- GaussRF(x = c(0,1, 1/n), y = c(0,1, 1/n), model = "stable",
                grid = TRUE, gridtriple = TRUE,
                param = c(mean=0, variance=1, nugget=0, scale=1, kappa=1))
par(mfrow=c(2,2))

# plots for 4 sliding windows (2 horizontal, 2 vertical)
fd2d <- fd.estimate(rf2d, methods="filter1",
        window.size = 100, step.size=100, plot.loglog = TRUE)

## End(Not run)

