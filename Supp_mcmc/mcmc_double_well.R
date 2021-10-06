
# mcmc_double_well.R
#
# Author: Bob Verity
# Date: 2021-10-06
#
# Purpose:
# Demonstrates how to implement PT in the mcmc package to solve the double-well
# problem.
#
# ------------------------------------------------------------------

# load mcmc package
library(mcmc)

# define model and MCMC parameters
gamma <- 20
burnin <- 1e3
samples <- 1e4
rungs <- 10

# define general log-joint function over state space of models and parameters
ludfun <- function(state, log.pseudo.prior) {
  
  # unpack thermodynamic power
  rung <- state[1]
  beta <- (rung - 1) / (rungs - 1)
  
  # unpack parameters
  params <- state[-1]
  mu <- params[1]
  
  # calculate power posterior
  loglike <- -gamma*(mu^2 - 1)^2
  logprior <- dunif(mu, min = -2, max = 2, log = TRUE)
  
  # calculate joint distribution in log space, with likelihood raised to the
  # power beta
  ret <- beta*loglike + logprior
  
  return(ret)
}

# make matrix of neighboring models, in our case between adjacent rungs
neighbors <- matrix(FALSE, rungs, rungs)
diag(neighbors[-1,]) <- diag(neighbors[,-1]) <- TRUE

# initialise state of chains
state_initial <- matrix(0, rungs, 1)

# run MCMC
out <- temper(obj = ludfun,
              initial = state_initial,
              neighbors = neighbors,
              nbatch = burnin + samples,
              log.pseudo.prior = 0,
              parallel = TRUE,
              scale = 0.3)

# get MC acceptance rates
mc_accept <- (diag(out$accepti[-1,]) + diag(out$accepti[,-1]))/2
plot(mc_accept)

# get draws from cold chain
cold_draws <- out$batch[-(1:burnin),rungs,]

# trace plot of draws from cold chain
plot(cold_draws)
