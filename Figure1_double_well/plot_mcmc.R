
# plot_mcmc.R
#
# Author: Bob Verity
# Date: 2019-11-29
#
# Purpose:
# Read in MCMC results from drjacoby and STAN for contrived double-well problem,
# and produce plots.
#
# ------------------------------------------------------------------

# load drjacoby from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/MCMC/drjacoby")

# load other packages
library(ggplot2)
library(cowplot)

# double well function
double_well <- function(x, gamma, return_log = FALSE) {
  ret <- -gamma*(x^2-1)^2
  if (!return_log) {
    ret <- exp(ret)
  }
  return(ret)
}

# double well function, normalised to integrate to 1
double_well_norm <- function(x, gamma, x_range, return_log = FALSE) {
  area <- integrate(double_well, lower = x_range[1], upper = x_range[2], gamma = gamma)$value
  ret <- double_well(x = x, gamma = gamma, return_log = return_log)
  if (return_log) {
    ret <- ret - log(area)
  } else {
    ret <- ret/area
  }
  return(ret)
}

# ------------------------------------------------------------------

# define plotting parameters
bw <- 0.05  # bin width
trace_size <- 0.3  # size of lines in trace plots
dens_size <- 0.7  # size of lines in density plot
dens_ylim <- c(0, 0.07)  # y-axis limits of density plot
col_vec <- drjacoby_cols()[c(2,4)]
col_dens <- "black"

# read in model parameters
params <- readRDS("Figure1_double_well/output/model_params.rds")
gamma_vec <- params$gamma_vec
mu_range <- params$mu_range
mu_vec <- seq(mu_range[1], mu_range[2], l = 1001)

# define double-well curves
df_true <- do.call(rbind, mapply(function(i) {
  ret <- data.frame(mu = mu_vec,
                    y = double_well_norm(mu_vec, gamma_vec[i], mu_range),
                    rep = sprintf("gamma = %s", gamma_vec[i]))
  ret$y <- 1 / (diff(mu_range)*bw) * ret$y / sum(ret$y)
  ret
}, 1:3, SIMPLIFY = FALSE))
df_true$rep <- factor(df_true$rep, levels = sprintf("gamma = %s", gamma_vec))

# read in MCMC output
drj_list <- readRDS("Figure1_double_well/output/drj_list.rds")
stan_list <- readRDS("Figure1_double_well/output/stan_list.rds")

# get STAN data into plotting dataframe
df_stan <- do.call(rbind, mapply(function(i) {
  ret <- stan_list[[i]]
  ret$method <- "STAN"
  ret$rep = sprintf("gamma = %s", gamma_vec[i])
  ret
}, 1:3, SIMPLIFY = FALSE))
df_stan$rep <- factor(df_stan$rep, levels = sprintf("gamma = %s", gamma_vec))

# get effective sample sizes
ess_stan <- mapply(function(x) round(coda::effectiveSize(x)), split(df_stan$mu, df_stan$rep))
ess_drj <- mapply(function(x) x$diagnostics$ess, drj_list)

# get drjacoby data into plotting dataframe
rungs <- max(drj_list[[1]]$output$rung)
df_drj <- do.call(rbind, mapply(function(i) {
  ret <- drj_list[[i]]$output
  ret <- subset(ret, phase == "sampling" & rung == rungs)
  ret$iteration <- seq_along(ret$iteration)
  ret$method <- "drj"
  ret$rep = sprintf("gamma = %s", gamma_vec[i])
  ret
}, 1:3, SIMPLIFY = FALSE))
df_drj$rep <- factor(df_drj$rep, levels = sprintf("gamma = %s", gamma_vec))

# ----------------------------

# produce base plot
#plot_base <- ggplot() + theme_bw() + theme(strip.background = element_blank(),
#                                           strip.text = element_blank())
plot_base <- ggplot() + theme_bw() + theme(strip.background = element_rect(fill = grey(1)))


# produce density plot
plot_dist <- plot_base +
  #geom_line(aes(x = mu, y = y), col = col_vec[1], size = dens_size, data = df_true) +
  geom_area(aes(x = mu, y = y), col = NA, fill = col_dens, size = dens_size, data = df_true) +
  stat_bin(aes(x = mu, y = (..count..)/sum(..count..)),
           geom = "step", size = dens_size, binwidth = bw, color = col_vec[1], data = df_stan) +
  stat_bin(aes(x = mu, y = (..count..)/sum(..count..)),
           geom = "step", size = dens_size, binwidth = bw, color = col_vec[2], data = df_drj) +
  scale_x_continuous(limits = mu_range, expand = c(0, 0)) +
  scale_y_continuous(limits = dens_ylim, expand = c(0, 0)) +
  facet_wrap(. ~rep, ncol = 1, strip.position = "left") +
  xlab("x") + ylab("probability") +
  ggtitle("posterior density")


# produce trace plot from STAN
plot_trace1 <- plot_base +
  geom_point(aes(x = iteration, y = mu), size = trace_size, color = col_vec[1], data = df_stan) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = mu_range, expand = c(0, 0)) +
  facet_wrap(. ~rep, ncol = 1, strip.position = "left") +
  xlab("iteration") + ylab("x") +
  ggtitle(expression(paste(italic("rstan"), " trace")))


# produce trace plot from drjacoby
plot_trace2 <- plot_base +
  geom_point(aes(x = iteration, y = mu), size = trace_size, color = col_vec[2], data = df_drj) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = mu_range, expand = c(0, 0)) +
  facet_wrap(. ~rep, ncol = 1, strip.position = "left") +
  xlab("iteration") + ylab("x") +
  ggtitle(expression(paste(italic("drjacoby"), " trace")))

# arrange in grid
plot_combined <- plot_grid(plot_trace1,  plot_trace2, plot_dist, ncol = 3)
plot_combined

# save to file
ggsave("Figure1_double_well/double_well.png", plot_combined, width = 8, height = 5)
