
# .R
#
# Author: Bob Verity
# Date: 2020-12-09
#
# Purpose:
# (this is an example header)
#
# ------------------------------------------------------------------
set.seed(9)

# define double well log-likelihood, up to constant of proportionality
f1 <- function(x, gamma) {
  -gamma*(x^2 - 1)^2
}

# define model and MCMC parameters
gamma <- 6
iterations <- 1e2
rungs <- 11
prop_sd <- 0.3
init_x <- -1
couping_on <- TRUE

# initialise
x <- runif(rungs, -3, 3)
x[1] <- init_x
beta_vec <- seq(1, 0, l = rungs)
ll <- f1(x, gamma)

# objects for storing results
cold_chain <- rep(NA, iterations)
track_chain <- 1
track_full <- matrix(NA, iterations, 2)
mc_accept <- rep(0, rungs - 1)

# MCMC
for (i in 1:iterations) {

  # store results
  cold_chain[i] <- x[1]
  track_full[i,] <- c(x[track_chain], track_chain)

  # basic MH
  x_prop <- x + rnorm(rungs, sd = prop_sd)
  w <- which((x_prop < -3) | (x_prop > 3))
  x_prop[w] <- x[w]
  ll_prop <- f1(x_prop, gamma)
  MH <- exp(beta_vec*(ll_prop - ll))
  w <- which(runif(rungs) < MH)
  if (any(w)) {
    x[w] <- x_prop[w]
    ll[w] <- ll_prop[w]
  }

  # coupling
  if (couping_on) {
    for (j in 2:rungs) {
      MC <- (beta_vec[j-1]*ll[j] + beta_vec[j]*ll[j-1]) - (beta_vec[j-1]*ll[j-1] + beta_vec[j]*ll[j])
      if (runif(1) < exp(MC)) {
        tmp <- x[j]
        x[j] <- x[j-1]
        x[j-1] <- tmp

        tmp <- ll[j]
        ll[j] <- ll[j-1]
        ll[j-1] <- tmp

        # update coupling acceptance
        mc_accept[j-1] <- mc_accept[j-1] + 1

        # update track chain
        if (j == track_chain) {
          track_chain <- j-1
        } else if (j-1 == track_chain) {
          track_chain <- j
        }

      }
    }
  }

}

# plot cold chain
plot(cold_chain)
hist(cold_chain, breaks = seq(-3, 3, 0.1), probability = TRUE)

# plot MC acceptance
plot(rev(mc_accept) / iterations, ylim = c(0,1))

# plot full track chain
plot(track_full, type = 'l', xlab = "x", ylab = "rung")

# subset to first cross-over
w <- which((track_full[,1] > 0.5) & (track_full[,2] == 1))[1]
track_full <- track_full[1:w,]
plot(track_full, type = 'l', xlab = "x", ylab = "rung")

# save results to file
ret <- list(gamma = gamma,
            rungs = rungs,
            track_full = track_full)
saveRDS(ret, "Figure_infographic/output/sim_results.rds")
