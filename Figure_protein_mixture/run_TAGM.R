# run_TAGM
#
# Author: Bob Verity
# Date: 2021-01-14
#
# Purpose:
# Read in the protein data and run TAGM mixture model MCMC for a range of K.
# MCMC objects are saved to file, although not back into this repository due to
# large file sizes.
#
# ------------------------------------------------------------------
# packages

# load drjacoby from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/MCMC/drjacoby")
library(drjacoby)

set.seed(1)


# ------------------------------------------------------------------
# FUNCTIONS

# source C++ likelihood and prior
Rcpp::sourceCpp("Figure_protein_mixture/source/TAGM_loglike.cpp")

# takes matrix as input, converts to list format for use within Rcpp code
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# density of multivariate t-distribution
dmvt <- function(x, mean_vec, scale_mat, deg_freedom, return_log = FALSE) {
  d <- length(mean_vec)
  ret <- lgamma((d + deg_freedom) / 2) - lgamma(deg_freedom / 2) - (d / 2)*log(pi * deg_freedom) - 0.5 * determinant(scale_mat, logarithm = TRUE)$modulus[1] -
    (deg_freedom + d) / 2 * log(1 + 1 / deg_freedom * ((x - mean_vec) %*% solve(scale_mat) %*% t(x - mean_vec))[1])
  if (!return_log) {
    ret <- exp(ret)
  }
  ret
}

# ------------------------------------------------------------------
# SETUP

# read in raw data
dat <- readRDS("Figure_protein_mixture/output/hyperLOPIT2015_processed.rds")

# perform PCA
p <- prcomp(na.omit(dat))
p$percent_var <- p$sdev^2 / sum(p$sdev^2) * 100
x <- p$x[,1]
y <- p$x[,2]

# what percent of variance is explained by first 2 components?
sum(p$percent_var[1:2])

# get global mean and covariance of data
mu_0 <- colMeans(p$x[,1:2])
S_0 <- cov(p$x[,1:2])

# define parameters of global t-distribution
t_mean <- t(mu_0)
t_scale <- S_0 / 2
t_df <- 4

# get probability of each observation from global t-distribution
t_prob <- apply(cbind(x, y), 1, dmvt, mean_vec = t_mean, scale_mat = t_scale, deg_freedom = t_df)

# define MCMC parameters
burnin <- 1e3
samples <- 1e4
chains <- 3

# ------------------------------------------------------------------
# MAIN K LOOP

for (K in 2:15) {
  #K <- 2
  
  # define free parameters
  df_params <- NULL
  for (i in 1:K) {
    df_params <- rbind(df_params,
                       define_params(name = sprintf("mu%s_x", i), min = -2, max = 2,
                                     name = sprintf("mu%s_y", i), min = -2, max = 2,
                                     name = sprintf("sigma%s_x", i), min = 0, max = 1,
                                     name = sprintf("sigma%s_y", i), min = 0, max = 1,
                                     name = sprintf("rho%s", i), min = -0.8, max = 0.8))
  }
  df_params <- rbind(df_params,
                     data.frame(name = sprintf("w_gamma%s", 1:K), min = 0, max = Inf),
                     data.frame(name = "eps", min = 0, max = 1))
  
  # define list of misc objects
  misc_list <- list(K = K,
                    mu_0 = mu_0,
                    lambda_0 = 0.01,
                    nu_0 = 4,
                    S_0 = matrix_to_rcpp(S_0),
                    beta = 1,
                    u = 2,
                    v = 999,
                    t_prob = t_prob)
  
  # run MCMC with single rung
  message(sprintf("\nBASIC, K = %s", K))
  mcmc_basic <- run_mcmc(data = list(x = x,
                                     y = y),
                         df_params = df_params,
                         misc = misc_list,
                         loglike = "loglike",
                         logprior = "logprior",
                         burnin = burnin,
                         samples = samples,
                         rungs = 1,
                         chains = chains)
  
  #plot_par(mcmc_basic)
  
  # run MCMC with multiple rungs
  message(sprintf("\nRUNGS, K = %s", K))
  mcmc_rungs <- run_mcmc(data = list(x = x,
                                     y = y),
                         df_params = df_params,
                         misc = misc_list,
                         loglike = "loglike",
                         logprior = "logprior",
                         burnin = burnin,
                         samples = samples,
                         beta_manual = seq(0, 1, 0.025),
                         alpha = 3.0,
                         chains = chains)
  
  #plot_rung_loglike(mcmc_rungs)
  #plot_mc_acceptance(mcmc_rungs)
  #plot_par(mcmc_rungs)
  
  # save both to file
  saveRDS(mcmc_basic, file = sprintf("/Users/rverity/Desktop/protein_mixture_output/TAGM_basic_%s", K))
  saveRDS(mcmc_rungs, file = sprintf("/Users/rverity/Desktop/protein_mixture_output/TAGM_rungs_%s", K))
  
}  # end K loop
