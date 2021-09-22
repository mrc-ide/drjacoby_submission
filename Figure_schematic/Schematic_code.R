# set random seed
set.seed(1)

# define true parameter values
mu_true <- 3
sigma_true <- 2

data_list <- list(x = rnorm(10, mean = mu_true, sd = sigma_true))

df_params <- define_params(name = "mu", min = -10, max = 10,
                           name = "sigma", min = 0, max = Inf)

r_loglike <- function(params, data, misc) {
  sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
}

r_logprior <- function(params, misc) {
  dunif(params["mu"], min = -10, max = 10, log = TRUE) +
    dlnorm(params["sigma"], meanlog = 0, sdlog = 1.0, log = TRUE)
}

