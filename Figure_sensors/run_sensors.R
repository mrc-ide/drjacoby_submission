
library(drjacoby)
library(dplyr)
library(ggplot2)
library(tidyr)
library(furrr)
library(patchwork)

set.seed(123987)

# Define parameter data.frame
s <- seq(-4, 4, length.out = 4)
init_x <- rep(s, 4)
init_y <- rep(s, each = 4)
params_mcmc <- define_params(name = "x1", min = -10, max = 10, init = init_x, 
                             name = "y1", min = -10, max = 10, init = init_y, 
                             name = "x2", min = -10, max = 10, init = init_x, 
                             name = "y2", min = -10, max = 10, init = init_y)
params_pt_mcmc <- define_params(name = "x1", min = -10, max = 10, init = init_x[1],
                                name = "y1", min = -10, max = 10, init = init_y[1],
                                name = "x2", min = -10, max = 10, init = init_x[1],
                                name = "y2", min = -10, max = 10, init = init_y[1])

Rcpp::sourceCpp("Figure_sensors/sensors.cpp")

data <- list(
  y = c(1.5, sqrt(1^2 + 0.5^2), 1.4)
)

sensor_output <- run_mcmc(data = data,
                          df_params = params_mcmc,
                          loglike = "loglike",
                          logprior = "logprior",
                          burnin = 1e3,
                          samples = 5e4,
                          chains = 16)

# Run the mcmc with tempering
sensor_output_tempered <- run_mcmc(data = data,
                                   df_params = params_pt_mcmc,
                                   loglike = "loglike",
                                   logprior = "logprior",
                                   burnin = 1e3,
                                   samples = 5e4,
                                   chains = 1,
                                   GTI_pow = 3,
                                   rungs = 16)

plot_mc_acceptance(sensor_output_tempered)
saveRDS(sensor_output, "Figure_sensors/output/sensor_output.RDS")
saveRDS(sensor_output_tempered, "Figure_sensors/output/sensor_output_tempered.RDS")

# True possible positions of sensor 1
sensor11 <- c(1.5, 1)
sensor12 <- c(0, -0.5)

td1 <- sensor_output$output %>%
  filter(phase == "sampling") %>%
  mutate(type = "mcmc")
td2 <- sensor_output_tempered$output %>%
  filter(phase == "sampling") %>%
  mutate(type = "ptmcmc")
td <- bind_rows(td1, td2)

ggplot(td, aes(x = x2, fill = type)) + 
  geom_density(alpha = 0.2)
ggplot(td, aes(x = y2, fill = type)) + 
  geom_density(alpha = 0.2)

## Grid search #################################################################

# Note this will require a decent amount of memory (can reduce resolution)

res = 0.05
grid <- bind_rows(
  expand_grid(
    x1 = seq(-1, 2, res),
    y1 = seq(-1, 2, res),
    x2 = seq(-2, 3.5, res),
    y2 = seq(-2.5, 3, res)
  )
)
nrow(grid) / 1e6

ll_lp <- function(x1, y1, x2, y2, data){
  p <- c(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
  loglike(p, data, list()) + logprior(p, list())
}

future::plan(multicore)
out <- unlist(future_pmap(grid, ll_lp, data = data))
grid$p <- exp(out)

pd <- grid %>%
  group_by(x2, y2) %>%
  summarise(p = mean(p))

ggplot(pd, aes(x = x2, y = y2, fill = p))+
  geom_tile()+
  scale_fill_viridis_c()

sampled_grid <- grid[sample(1:nrow(grid), 100000, replace = TRUE, prob = grid$p),]
saveRDS(sampled_grid, "Figure_sensors/output/sampled_grid.RDS")
################################################################################

### Plotting ###################################################################
chain_order <- paste(c(13:16, 9:12, 5:8, 1:4))

chain_plot_data <- sensor_output$output %>%
  filter(phase == "sampling") %>%
  mutate(type = "MCMC") %>%
  mutate(chain = factor(chain, levels = chain_order))
starting_point_data <- data.frame(x2 = init_x, y2 = init_y, chain = 1:16) %>%
  mutate(chain = factor(chain, levels = chain_order))

chain_plot <- ggplot() +
  geom_point(data = slice_sample(sampled_grid, n = 2000), aes(x = x2, y = y2), col = "grey", size = 0.1) +
  geom_hex(data = chain_plot_data, aes(x = x2, y = y2, fill = ..density..), bins = 100) +
  geom_point(data = starting_point_data, aes(x = x2, y = y2), col = "darkcyan", fill = "cyan2", shape = 21, size = 2) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~ chain, nrow = 5) + 
  xlab("x") +
  ylab("y") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

pt_mcmc_plot_data <- sensor_output_tempered$output %>%
  filter(phase == "sampling") %>%
  mutate(type = "PT MCMC") %>%
  mutate(chain = factor(chain, levels = chain_order))
posterior_plot_data <- chain_plot_data %>%
  slice_sample(n = nrow(pt_mcmc_plot_data)) %>%
  bind_rows(pt_mcmc_plot_data)
posterior_plot_data <- posterior_plot_data[sample(nrow(posterior_plot_data)),]
# XY overlay
dxy <- ggplot() +
  geom_point(data = posterior_plot_data, aes(x = x2, y = y2, col = type), size = 0.5, alpha = 0.2) +
  geom_line(aes(x = 0, y = 0, linetype = "True", group = 1)) +
  scale_linetype_manual(values =  2) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.title = element_blank())

# X density
dx <- ggplot() + 
  geom_density(data = posterior_plot_data, aes(x = x2, fill = type, col = type), alpha = 0.3, bw = 0.1) +
  geom_density(data = sampled_grid, aes(x = x2), col = "black", lty = 2, bw = 0.1) +
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
  geom_density(data = posterior_plot_data, aes(x = y2, fill = type, col = type), alpha = 0.3, bw = 0.1) +
  geom_density(data = sampled_grid, aes(x = y2), col = "black", lty = 2, bw = 0.1) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

combo <- dxy + dy + dx + guide_area() + plot_layout(guides = 'collect')
sensor_plot <- chain_plot | combo

ggsave("Figure_sensors/sensors.png", sensor_plot, dpi = 600, height = 5, width = 9)
