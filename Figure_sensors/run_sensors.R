# Author: Pete Winskill
# Date: 2021-09027
#
# Purpose:
# Run MCMC in drjacoby for the sensor location example
#
# Inspired by the following work:
# Hyungsuk Tak, Xiao-Li Meng & David A. van Dyk (2018) A Repelling–Attracting
# Metropolis Algorithm for Multimodality, Journal of Computational and Graphical Statistics, 27:3,
# 479-490, DOI: 10.1080/10618600.2017.1415911
# ------------------------------------------------------------------


# Load packages
library(drjacoby)
library(dplyr)
library(ggplot2)
library(tidyr)
library(furrr)
library(patchwork)

# Set seed to reproduce results exactly
set.seed(12345678910)

# Define starting position grid
s <- seq(-4, 4, length.out = 4)
init_x <- rep(s, 4)
init_y <- rep(s, each = 4)

# Define parameters data.frame
params_mcmc <- define_params(name = "x3", min = -10, max = 10, init = init_x, 
                             name = "y3", min = -10, max = 10, init = init_y, 
                             name = "x4", min = -10, max = 10, init = init_x, 
                             name = "y4", min = -10, max = 10, init = init_y)
params_pt_mcmc <- define_params(name = "x3", min = -10, max = 10, init = init_x[1],
                                name = "y3", min = -10, max = 10, init = init_y[1],
                                name = "x4", min = -10, max = 10, init = init_x[1],
                                name = "y4", min = -10, max = 10, init = init_y[1])

# Source c++ log likelihood and log prior functions
Rcpp::sourceCpp("Figure_sensors/sensors.cpp")

# Define data. Observed distance between sensor 1 and 3, 2 and 3 and 3 and 4 
  # Sensors 1 and 2 have known locations at c(0, 1) and c(1, 0)
  # Based observed distance we know that sensor 3 can only be in two locations (1.5, 1) or (0, -0.5) 
  # with sensor 4 having possible rings of locations surrounding those coordinates
data <- list(
  y = c(1.5, sqrt(1^2 + 0.5^2), 1.4)
)

# Run standard MCMC
sensor_output <- run_mcmc(data = data,
                          df_params = params_mcmc,
                          loglike = "loglike",
                          logprior = "logprior",
                          burnin = 1e4,
                          samples = 1e5,
                          chains = 16)

# Run the mcmc with tempering
sensor_output_tempered <- run_mcmc(data = data,
                                   df_params = params_pt_mcmc,
                                   loglike = "loglike",
                                   logprior = "logprior",
                                   burnin = 1e4,
                                   samples = 1e5,
                                   chains = 1,
                                   alpha = 4,
                                   rungs = 16)
# View mc acceptance rates
plot_mc_acceptance(sensor_output_tempered)

# Save output
saveRDS(sensor_output, "Figure_sensors/output/sensor_output.RDS")
saveRDS(sensor_output_tempered, "Figure_sensors/output/sensor_output_tempered.RDS")

## Grid search #################################################################

# Note: this is memory hungry
res = 0.05
grid <- bind_rows(
  expand_grid(
    x3 = seq(-1, 2, res),
    y3 = seq(-1, 2, res),
    x4 = seq(-2, 3.5, res),
    y4 = seq(-2.5, 3, res)
  )
)

# A function to evaluate the log-likelihood + logprior
ll_lp <- function(x3, y3, x4, y4, data){
  p <- c(x3 = x3, y3 = y3, x4 = x4, y4 = y4)
  loglike(p, data, list()) + logprior(p, list())
}

# Apply to grid (parallelised)
future::plan(multicore)
out <- unlist(future_pmap(grid, ll_lp, data = data))
# Exponentiate:
grid$p <- exp(out)

# Sample the grid based on the probability
sampled_grid <- grid[sample(1:nrow(grid), 100000, replace = TRUE, prob = grid$p),]
saveRDS(sampled_grid, "Figure_sensors/output/sampled_grid.RDS")
################################################################################

### Plotting ###################################################################

sensor_output <- readRDS("Figure_sensors/output/sensor_output.RDS")
sensor_output_tempered <- readRDS("Figure_sensors/output/sensor_output_tempered.RDS")
sampled_grid <- readRDS("Figure_sensors/output/sampled_grid.RDS")

# Plotting individual chains
chain_order <- paste(c(13:16, 9:12, 5:8, 1:4))

chain_plot_data <- sensor_output$output %>%
  filter(phase == "sampling") %>%
  mutate(type = "Standard") %>%
  mutate(chain = factor(chain, levels = chain_order))
starting_point_data <- data.frame(x4 = init_x, y4 = init_y, chain = 1:16) %>%
  mutate(chain = factor(chain, levels = chain_order))

chain_plot <- ggplot() +
  geom_point(data = slice_sample(sampled_grid, n = 2000), aes(x = x4, y = y4), col = "grey", size = 0.1) +
  geom_hex(data = chain_plot_data, aes(x = x4, y = y4, fill = ..density..), bins = 100) +
  geom_point(data = starting_point_data, aes(x = x4, y = y4), col = "#f9c22e", fill = "#f9c22e", shape = 21, size = 2) +
  theme_bw() +
  scale_fill_viridis_c(option = "A", name = "Density") +
  facet_wrap(~ chain, nrow = 5) + 
  xlab("x") +
  ylab("y") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(size = 9), 
        legend.text = element_text(size = 7), 
        aspect.ratio = 1)

pt_mcmc_plot_data <- sensor_output_tempered$output %>%
  filter(phase == "sampling") %>%
  mutate(type = "PT") %>%
  mutate(chain = factor(chain, levels = chain_order))
posterior_plot_data <- chain_plot_data %>%
  slice_sample(n = nrow(pt_mcmc_plot_data)) %>%
  bind_rows(pt_mcmc_plot_data)
posterior_plot_data <- posterior_plot_data[sample(nrow(posterior_plot_data)),]
# XY overlay
dxy <- ggplot() +
  geom_point(data = slice_sample(posterior_plot_data, n = 50000), aes(x = x4, y = y4, col = type), size = 0.5, alpha = 0.2) +
  scale_colour_manual(values = c("#53b3cb", "#cc282b")) +
  scale_linetype_manual(values =  2) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
  xlab("x") +
  ylab("y") +
  theme(legend.title = element_blank())

# X density
dx <- ggplot() + 
  geom_density(data = sampled_grid, aes(x = x4), fill = "grey", col = NA, bw = 0.1) +
  geom_density(data = posterior_plot_data, aes(x = x4, col = type), alpha = 0.3, bw = 0.1, size = 1) +
  scale_colour_manual(values = c("#53b3cb", "#cc282b")) +
  theme_bw() +
  scale_y_reverse() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
# Y density
dy <- ggplot() + 
  geom_density(data = sampled_grid, aes(x = x4), fill = "grey", col = NA, bw = 0.1) +
  geom_density(data = posterior_plot_data, aes(x = x4, col = type), alpha = 0.3, bw = 0.1, size = 1) +
  scale_colour_manual(values = c("#53b3cb", "#cc282b")) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# Combine plots
combo <- dxy + dy + dx + guide_area() + plot_layout(guides = 'collect')
sensor_plot <- (chain_plot | combo) + plot_annotation(tag_levels = "A") + plot_layout(widths = c(0.7, 0.7))

# Save plot
ggsave("Figure_sensors/sensors.png", sensor_plot, dpi = 600, height = 4, width = 8)


# Sensor 3 plot for SI
sensor_3 <- plot_cor(sensor_output, "x3", "y3") +
  xlab("x") + 
  ylab("y")
ggsave("Figure_sensors/SI_sensors_3.png", sensor_3, dpi = 600, height = 4, width = 4)
