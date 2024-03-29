
# run_mcmc.R
#
# Author: Bob Verity
# Date: 2020-12-08
#
# Purpose:
# Compare sampling efficiency of various MCMC methods by a series of metrics,
# including run-times and effective sample sizes (ESS).
#
# ------------------------------------------------------------------

# load drjacoby from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/MCMC/drjacoby")
library(drjacoby)

# load other packages
library(rstan)
library(magrittr)

# double well function
double_well <- function(x, gamma, return_log = FALSE) {
  ret <- -gamma*(x^2 - 1)^2
  if (!return_log) {
    ret <- exp(ret)
  }
  return(ret)
}

# ------------------------------------------------------------------

set.seed(2)

# define model parameters
gamma_vec <- c(1, 5, 10, 15, 20)
mu_range <- c(-2, 2)

# define MCMC parameters
burnin <- 1e4
samples <- 1e5
rungs <- 10

# define parameters dataframe
df_params <- define_params(name = "mu", min = mu_range[1], max = mu_range[2])

# source C++ likelihood and prior
Rcpp::sourceCpp("Table_double_well/source/doublewell_loglike_logprior.cpp")

# define object for storing results
res_df <- data.frame(gamma = gamma_vec,
                     stan_warmup = NA,
                     stan_sample = NA,
                     stan_eval = NA,
                     stan_ESS = NA,
                     stan_mean = NA,
                     stan_Q2.5 = NA,
                     stan_Q97.5 = NA,
                     stan_long_warmup = NA,
                     stan_long_sample = NA,
                     stan_long_eval = NA,
                     stan_long_ESS = NA,
                     stan_long_mean = NA,
                     stan_long_Q2.5 = NA,
                     stan_long_Q97.5 = NA,
                     naive_internal = NA,
                     naive_eval = NA,
                     naive_ESS = NA,
                     naive_mean = NA,
                     naive_Q2.5 = NA,
                     naive_Q97.5 = NA,
                     naive_long_internal = NA,
                     naive_long_eval = NA,
                     naive_long_ESS = NA,
                     naive_long_mean = NA,
                     naive_long_Q2.5 = NA,
                     naive_long_Q97.5 = NA,
                     PT_internal = NA,
                     PT_eval = NA,
                     PT_ESS = NA,
                     PT_mean = NA,
                     PT_Q2.5 = NA,
                     PT_Q97.5 = NA)

# loop through values of gamma
for (i in seq_along(gamma_vec)) {
  message(sprintf("gamma = %s", gamma_vec[i]))
  
  # ----------------------------
  # STAN
  
  # run MCMC
  t0 <- Sys.time()
  fit <- stan(file = 'Figure_double_well/stan_models/double_well.stan',
              data = list(gamma = gamma_vec[i]),
              chains = 1,
              seed = 1,
              warmup = burnin,
              iter = burnin + samples)
  t1 <- Sys.time() - t0
  
  # get various run-times
  res_df$stan_warmup[i] <- get_elapsed_time(fit)[1]
  res_df$stan_sample[i] <- get_elapsed_time(fit)[2]
  res_df$stan_eval[i] <- t1
  
  # get ESS and postior estimates
  fit_summary <- as.data.frame(summary(fit)$summary)
  res_df$stan_ESS[i] <- fit_summary$n_eff[1]
  res_df$stan_mean[i] <- fit_summary$mean[1]
  res_df$stan_Q2.5[i] <- fit_summary$`2.5%`[1]
  res_df$stan_Q97.5[i] <- fit_summary$`97.5%`[1]
  
  # ----------------------------
  # STAN for longer
  
  # run MCMC
  t0 <- Sys.time()
  fit_long <- stan(file = 'Figure_double_well/stan_models/double_well.stan',
                   data = list(gamma = gamma_vec[i]),
                   chains = 1,
                   seed = 1,
                   warmup = burnin * rungs,
                   iter = (burnin + samples) * rungs)
  t1 <- Sys.time() - t0
  
  # get various run-times
  res_df$stan_long_warmup[i] <- get_elapsed_time(fit_long)[1]
  res_df$stan_long_sample[i] <- get_elapsed_time(fit_long)[2]
  res_df$stan_long_eval[i] <- t1
  
  # get ESS and posterior estimates
  fit_summary <- as.data.frame(summary(fit_long)$summary)
  res_df$stan_long_ESS[i] <- fit_summary$n_eff[1]
  res_df$stan_long_mean[i] <- fit_summary$mean[1]
  res_df$stan_long_Q2.5[i] <- fit_summary$`2.5%`[1]
  res_df$stan_long_Q97.5[i] <- fit_summary$`97.5%`[1]
  
  # ----------------------------
  # Naive (drjacoby without rungs)
  
  # run MCMC
  t0 <- Sys.time()
  drj_simple <- run_mcmc(data = list(gamma = gamma_vec[i]),
                         df_params = df_params,
                         loglike = "loglike",
                         logprior = "logprior",
                         burnin = burnin,
                         samples = samples,
                         chains = 1,
                         rungs = 1)
  t1 <- Sys.time() - t0
  
  # get various run-times
  res_df$naive_internal[i] <- sum(drj_simple$diagnostics$run_time$seconds)
  res_df$naive_eval[i] <- t1
  
  # get ESS and posterior estimates
  res_df$naive_ESS[i] <- drj_simple$diagnostics$ess[1]
  m <- drj_simple$output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::pull("mu") %>%
    mean()
  q <- drj_simple$output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::pull("mu") %>%
    quantile(probs = c(0.025, 0.975))
  res_df$naive_mean[i] <- m
  res_df$naive_Q2.5[i] <- q[1]
  res_df$naive_Q97.5[i] <- q[2]
  
  # ----------------------------
  # Naive for longer
  
  # run MCMC
  t0 <- Sys.time()
  drj_simple_long <- run_mcmc(data = list(gamma = gamma_vec[i]),
                              df_params = df_params,
                              loglike = "loglike",
                              logprior = "logprior",
                              burnin = burnin * rungs,
                              samples = samples *rungs,
                              chains = 1,
                              rungs = 1)
  t1 <- Sys.time() - t0
  
  # get various run-times
  res_df$naive_long_internal[i] <- sum(drj_simple_long$diagnostics$run_time$seconds)
  res_df$naive_long_eval[i] <- t1
  
  # get ESS and posterior estimates
  res_df$naive_long_ESS[i] <- drj_simple_long$diagnostics$ess[1]
  m <- drj_simple_long$output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::pull("mu") %>%
    mean()
  q <- drj_simple_long$output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::pull("mu") %>%
    quantile(probs = c(0.025, 0.975))
  res_df$naive_long_mean[i] <- m
  res_df$naive_long_Q2.5[i] <- q[1]
  res_df$naive_long_Q97.5[i] <- q[2]
  
  # ----------------------------
  # Parallel Tempering (drjacoby with rungs)
  
  # run MCMC
  t0 <- Sys.time()
  drj_mcmc_rungs <- run_mcmc(data = list(gamma = gamma_vec[i]),
                             df_params = df_params,
                             loglike = "loglike",
                             logprior = "logprior",
                             burnin = burnin,
                             samples = samples,
                             chains = 1,
                             rungs = rungs)
  t1 <- Sys.time() - t0
  
  # get various run-times
  res_df$PT_internal[i] <- sum(drj_mcmc_rungs$diagnostics$run_time$seconds)
  res_df$PT_eval[i] <- t1
  
  # get ESS and posterior estimates
  res_df$PT_ESS[i] <- drj_mcmc_rungs$diagnostics$ess[1]
  m <- drj_mcmc_rungs$output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::pull("mu") %>%
    mean()
  q <- drj_mcmc_rungs$output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::pull("mu") %>%
    quantile(probs = c(0.025, 0.975))
  res_df$PT_mean[i] <- m
  res_df$PT_Q2.5[i] <- q[1]
  res_df$PT_Q97.5[i] <- q[2]
  
}

#View(res_df)

# create simplified version of full results
res_simple <- res_df %>%
  dplyr::mutate(stan_long_internal = stan_long_warmup + stan_long_sample) %>%
  dplyr::select(gamma,
                naive_long_internal, naive_long_eval, naive_long_ESS, naive_long_mean, naive_long_Q2.5, naive_long_Q97.5,
                stan_long_internal, stan_long_eval, stan_long_ESS, stan_long_mean, stan_long_Q2.5, stan_long_Q97.5,
                PT_internal, PT_eval, PT_ESS, PT_mean, PT_Q2.5, PT_Q97.5) %>%
  signif(digits = 3)

# write to file
write.csv(res_df, "Table_double_well/output/table_double_well.csv",
          quote = FALSE, row.names = FALSE)
write.csv(res_simple, "Table_double_well/output/table_double_well_simple.csv",
          quote = FALSE, row.names = FALSE)



