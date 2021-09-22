# run_mcmc.R
#
# Author: Pete Winskill
# Date: 2020-12-08
#
# Purpose:
# Run MCMC in drjacoby for the sensor location example
#
# References:
# Hyungsuk Tak, Xiao-Li Meng & David A. van Dyk (2018) A Repelling–Attracting
# Metropolis Algorithm for Multimodality, Journal of Computational and Graphical Statistics, 27:3,
# 479-490, DOI: 10.1080/10618600.2017.1415911
# ------------------------------------------------------------------

# Load packages
library(drjacoby)

set.seed(98765)

# Distance between sensor subscripts
location = data.frame(i = c(2, 1, 2, 1, 3, 1, 3),
                      j = c(3, 4, 4, 5, 5, 6, 6))

# Observed distances (the order matters for use in the likelihood)
data <- list(y = c(0.2970, 0.9266, 0.8524, 0.6103, 0.2995, 0.3631, 0.5656))

# Define parameter data.frame
params1 <- define_params(name = "x1", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "y1", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "x2", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "y2", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "x3", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "y3", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "x4", min = -10, max = 10, init = runif(20, -1, 1),
                        name = "y4", min = -10, max = 10, init = runif(20, -1, 1))
# duplicate parameter data.frame with first init values for PT
params2 <- params1
params2$init <- sapply(params2$init, function(x){x[1]})

# Pre-calculate grid of unique sensor pairs
d <- expand.grid(i = 1:6, j = 1:6)
d <- d[d$i != d$j,]
d <- d[d$j > d$i,]
d$known <- c(0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0)
pairs <- list(d = d)

# log-Likelihood
ll <- function(params, data, misc){
  
  # Coordinates of sensors (4 unknown that we estimate and 2 known locations)
  coors <- matrix(c(params, 0.5, 0.3, 0.3 ,0.7), ncol = 2, byrow = TRUE)
  # Pairs
  d <- misc$d
  
  # Predicted Euclidean distances between pairs of sensors
  pred_y <- as.matrix(dist(coors))[cbind(d$i, d$j)]
  
  # Comparing observed and estimated distances between sensors
  ll1 <- sum(dnorm(data$y, pred_y[d$known == 1], 0.02, log = TRUE))
  
  # Assume the the probability of observing a nearby sensor falls exponentially with distance
  prob_y <- exp(-((pred_y^2) / (2 * 0.3 ^ 2)))
  ll2 <- sum(dbinom(d$known, 1, prob_y, log = TRUE))
  
  ll1 + ll2
}

# Simple diffuse prior on each coordinate
lp <- function(params, misc){
  sum(dnorm(params, 0, 10, log = TRUE))
}


# Run the mcmc with tempering
sensor_output <- run_mcmc(data = data,
                          df_params = params1,
                          loglike = ll,
                          logprior = lp,
                          misc = pairs,
                          burnin = 1e3,
                          samples = 1e4,
                          chains = 20)

# Run the mcmc with tempering
sensor_output_tempered <- run_mcmc(data = data,
                                   df_params = params2,
                                   loglike = ll,
                                   logprior = lp,
                                   misc = pairs,
                                   burnin = 1e3,
                                   samples = 1e4,
                                   chains = 1,
                                   GTI_pow = 5,
                                   rungs = 30)

# Save output
saveRDS(sensor_output, "Figure_sensors/output/sensor_output.RDS")
saveRDS(sensor_output_tempered, "Figure_sensors/output/sensor_output_tempered.RDS")
