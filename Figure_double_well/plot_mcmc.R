
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
library(RColorBrewer)

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
bw <- c(0.05, 0.02, 0.01)  # smoothing bandwidth of density plot
dens_size <- 0.7           # size of lines in density plot
dens_ylim <- c(0, 0.07)    # y-axis limits of density plot
col_vec <- c("#cc282b", "#344cd5", "#53b3cb")
col_dens <- grey(0.6)

# read in model parameters
params <- readRDS("Figure_double_well/output/model_params.rds")
gamma_vec <- params$gamma_vec
mu_range <- params$mu_range
mu_vec <- seq(mu_range[1], mu_range[2], l = 1001)

# define double-well curves
model_levels <- sprintf("gamma = %s", gamma_vec)
df_true <- do.call(rbind, mapply(function(i) {
  ret <- data.frame(mu = mu_vec,
                    y = double_well_norm(mu_vec, gamma_vec[i], mu_range),
                    model = model_levels[i])
  ret$y <- ret$y / sum(ret$y)
  ret
}, 1:3, SIMPLIFY = FALSE))
df_true$model <- factor(df_true$model, levels = model_levels)

# read in MCMC output
drj_list_simple <- readRDS("Figure_double_well/output/drj_list_simple.rds")
drj_list_rungs <- readRDS("Figure_double_well/output/drj_list_rungs.rds")
stan_list <- readRDS("Figure_double_well/output/stan_list.rds")

# ----------------------------

# get STAN data into plotting dataframe
df_stan <- do.call(rbind, mapply(function(i) {
  mu <- stan_list[[i]]$mu
  d <- density(mu, from = -2, to = 2, bw = bw[i])
  data.frame(mu = d$x, y = d$y, model = model_levels[i])
}, 1:3, SIMPLIFY = FALSE))
df_stan$model <- factor(df_stan$model, levels = model_levels)

# get simple drjacoby data into plotting dataframe
df_drj_simple <- do.call(rbind, mapply(function(i) {
  x <- subset(drj_list_simple[[i]]$output, phase == "sampling")
  mu <- x$mu
  d <- density(mu, from = -2, to = 2, bw = bw[i])
  data.frame(mu = d$x, y = d$y, model = model_levels[i])
}, 1:3, SIMPLIFY = FALSE))
df_drj_simple$model <- factor(df_drj_simple$model, levels = model_levels)

# get tempered drjacoby data into plotting dataframe
df_drj_rungs <- do.call(rbind, mapply(function(i) {
  x <- subset(drj_list_rungs[[i]]$output, phase == "sampling")
  mu <- x$mu
  d <- density(mu, from = -2, to = 2, bw = bw[i])
  data.frame(mu = d$x, y = d$y, model = model_levels[i])
}, 1:3, SIMPLIFY = FALSE))
df_drj_rungs$model <- factor(df_drj_rungs$model, levels = model_levels)

# produce density plot
plot_dist <- ggplot() + theme_bw() +
  geom_area(aes(x = mu, y = 250*y), col = NA, fill = col_dens, size = dens_size, data = df_true) +
  geom_line(aes(x = mu, y = y), col = col_vec[1], data = df_drj_simple) +
  geom_line(aes(x = mu, y = y), col = col_vec[2], data = df_stan) +
  geom_line(aes(x = mu, y = y), col = col_vec[3], data = df_drj_rungs) +
  scale_x_continuous(limits = mu_range, expand = c(0, 0)) +
  facet_wrap(. ~model, ncol = 1, strip.position = "left", scales = "free_y") +
  xlab("mu") + ylab("Probability") +
  ggtitle("Posterior density") +
  theme(strip.background = element_rect(fill = "white"),
        panel.spacing.y = unit(1, "lines")) +
  labs(tag = "B")

# ----------------------------

# get posterior draws into single dataframe over models
df_trace <- do.call(rbind, mapply(function(i) {
  
  # extract results for this gamma
  df_standard <- drj_list_simple[[i]]$output %>%
    dplyr::mutate(model = "Basic MH")
  df_PT <- drj_list_rungs[[i]]$output %>%
    dplyr::mutate(model = "PT")
  df_stan <- stan_list[[i]] %>%
    dplyr::mutate(model = "HMC")
  
  # combine
  df_standard %>%
    dplyr::bind_rows(df_PT) %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::select(iteration, mu, model) %>%
    dplyr::mutate(iteration = iteration - min(iteration) + 1) %>%
    dplyr::bind_rows(df_stan) %>%
    dplyr::mutate(gamma = sprintf("gamma = %s", gamma_vec[i]),
                  model = factor(model, levels = c("Basic MH", "HMC", "PT")))
}, 1:3, SIMPLIFY = FALSE))

# produce combined trace plot
plot_trace <- df_trace %>%
  dplyr::mutate(gamma = factor(gamma, levels = sprintf("gamma = %s", gamma_vec))) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = mu, col = model), size = 0.1) +
  facet_grid(gamma ~ model, switch = "y") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = mu_range, expand = c(0, 0)) +
  scale_color_manual(values = col_vec) +
  xlab("Iteration") +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.spacing.y = unit(1, "lines")) +
  ggtitle("Trace plots") +
  labs(tag = "A")

# ----------------------------

# produce combined plot
plot_combined <- plot_grid(plot_trace, plot_dist, nrow = 1, rel_widths = c(3, 1))
plot_combined

# save to file
ggsave("Figure_double_well/double_well.png", plot_combined, width = 10, height = 5)
