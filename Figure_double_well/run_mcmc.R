
# run_mcmc.R
#
# Author: Bob Verity
# Date: 2020-12-08
#
# Purpose:
# Run MCMC in drjacoby and STAN under double-well model.
#
# ------------------------------------------------------------------

# load drjacoby from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/MCMC/drjacoby")

# load other packages
library(rstan)

# double well function
double_well <- function(x, gamma, return_log = FALSE) {
  ret <- -gamma*(x^2-1)^2
  if (!return_log) {
    ret <- exp(ret)
  }
  return(ret)
}

# ------------------------------------------------------------------

set.seed(2)

# define model parameters
gamma_vec <- c(1, 5, 10)
mu_range <- c(-2, 2)

# save parameters to file
saveRDS(list(gamma_vec = gamma_vec,
             mu_range = mu_range),
        "Figure_double_well/output/model_params.rds")

# define R logprior function
r_logprior <- function(params, misc) {
  dunif(params["mu"], mu_range[1], mu_range[2], log = TRUE)
}

# define MCMC parameters
burnin <- 1e3
samples <- 1e4
rungs <- 10

# define parameters dataframe
df_params <- define_params(name = "mu", min = -2, max = 2)

# loop through values of gamma
drj_list_simple <- drj_list_rungs <- stan_list <- list()
for (i in 1:3) {
  
  # run STAN MCMC
  fit <- stan(file = 'Figure_double_well/stan_models/double_well.stan',
              data = list(gamma = gamma_vec[i]),
              chains = 1,
              seed = 1,
              warmup = burnin,
              iter = burnin + samples)
  fit_extract <- rstan::extract(fit)
  stan_mcmc <- data.frame(iteration = 1:length(fit_extract$mu),
                          mu = fit_extract$mu)
  
  
  # define R loglike function
  r_loglike <- function(theta, data, misc) {
    double_well(theta[1], gamma = gamma_vec[i], return_log = TRUE)
  }
  
  # run drjacoby MCMC without rungs
  drj_mcmc_simple <- run_mcmc(data = list(x = -1),
                              df_params = df_params,
                              loglike = r_loglike,
                              logprior = r_logprior,
                              burnin = burnin,
                              samples = samples,
                              chains = 1,
                              rungs = 1)
  
  # run drjacoby MCMC with rungs
  drj_mcmc_rungs <- run_mcmc(data = list(x = -1),
                             df_params = df_params,
                             loglike = r_loglike,
                             logprior = r_logprior,
                             burnin = burnin,
                             samples = samples,
                             chains = 1,
                             rungs = 10)
  
  #plot_mc_acceptance(drj_mcmc_rungs)
  
  # save to list
  stan_list[[i]] <- stan_mcmc
  drj_list_simple[[i]] <- drj_mcmc_simple
  drj_list_rungs[[i]] <- drj_mcmc_rungs
  
}

# save to file
saveRDS(stan_list, "Figure_double_well/output/stan_list.rds")
saveRDS(drj_list_simple, "Figure_double_well/output/drj_list_simple.rds")
saveRDS(drj_list_rungs, "Figure_double_well/output/drj_list_rungs.rds")

