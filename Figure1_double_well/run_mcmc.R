
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

set.seed(1)

# define model parameters
gamma_vec <- c(1, 6.25, 10)
mu_range <- c(-2, 2)

# save parameters to file
saveRDS(list(gamma_vec = gamma_vec,
             mu_range = mu_range),
        "Figure1_double_well/output/model_params.rds")

# define R logprior function
r_logprior <- function(params, misc) {
  dunif(params["mu"], mu_range[1], mu_range[2], log = TRUE)
}

# define MCMC parameters
burnin <- 1e3
samples <- 1e3
chains <- 1
rungs <- 1

# define parameters dataframe
df_params <- define_params(name = "mu", min = -2, max = 2)

# loop through values of gamma
drj_list <- stan_list <- list()
for (i in 1:3) {

  # run STAN MCMC
  fit <- stan(file = 'Figure1_double_well/stan_models/double_well.stan',
              data = list(gamma = gamma_vec[i]), chains = 1)
  fit_extract <- extract(fit)
  stan_mcmc <- data.frame(iteration = 1:length(fit_extract$mu),
                          mu = fit_extract$mu)


  # define R loglike function
  r_loglike <- function(theta, data, misc) {
    double_well(theta[1], gamma = gamma_vec[i], return_log = TRUE)
  }

  # run drjacoby MCMC
  drj_mcmc <- run_mcmc(data = list(x = -1),
                       df_params = df_params,
                       loglike = r_loglike,
                       logprior = r_logprior,
                       burnin = burnin,
                       samples = samples,
                       chains = chains,
                       rungs = 20)

  # save to list
  stan_list[[i]] <- stan_mcmc
  drj_list[[i]] <- drj_mcmc

}

# save to file
saveRDS(stan_list, "Figure1_double_well/output/stan_list.rds")
saveRDS(drj_list, "Figure1_double_well/output/drj_list.rds")

