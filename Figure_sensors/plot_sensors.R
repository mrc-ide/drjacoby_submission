# plot_sensors.R
#
# Author: Pete Winskill
# Date: 2020-12-08
#
# Purpose:
# Plotting MCMC in drjacoby for the sensor location example
#
# Reference:
# Hyungsuk Tak, Xiao-Li Meng & David A. van Dyk (2018) A Repelling–Attracting
# Metropolis Algorithm for Multimodality, Journal of Computational and Graphical Statistics, 27:3,
# 479-490, DOI: 10.1080/10618600.2017.1415911
# ------------------------------------------------------------------


# Load packages
library(drjacoby)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load mcmc output
sensor_output <- readRDS("Figure_sensors/output/sensor_output.RDS")
sensor_output_tempered <- readRDS("Figure_sensors/output/sensor_output_tempered.RDS")
sampled_grid.RDS <- readRDS("Figure_sensors/output/sampled_grid.RDS")

# Which sensor to plot:
sensor <- 2
xname <- paste0("x", sensor)
yname <- paste0("y", sensor)

chain_order <- paste(c(21, 22, 23, 24, 25, 16, 17, 18, 19, 20, 11, 12, 13, 14, 15, 6, 7, 8, 9, 10, 1, 2, 3, 4,5))

# Extract starting positions
sp <- data.frame(x = unlist(params1[params1$name == xname, "init"]),
                 y = unlist(params1[params1$name == yname, "init"])) %>%
  mutate(chain = 1:25,
         chain = factor(chain, levels = chain_order))

mcmc <- sensor_output$output %>%
  mutate(type = "MCMC") %>%
  filter(phase == "sampling") %>%
  select(all_of(c("chain", "type", xname, yname))) %>%
  mutate(chain = factor(chain, levels = chain_order))
colnames(mcmc) <- c("chain", "type", "x", "y")

pt <- sensor_output_tempered$output %>%
  mutate(type = "PT MCMC") %>%
  filter(phase == "sampling") %>%
  select(all_of(c("chain", "type", xname, yname))) %>%
  mutate(chain = factor(chain, levels = chain_order))
colnames(pt) <- c("chain", "type", "x", "y")

# Make sure we plot the same number of samples
combined <- bind_rows(mcmc, pt) %>%
  select(-chain)
# Mix up samples for plotting
combined_sub <- bind_rows(slice_sample(mcmc, n = nrow(pt)), pt) %>%
  select(-chain)
combined_sub <- combined_sub[sample(nrow(combined_sub)),]


# Plot grid of MCMC chains
dgrid <- ggplot() + 
  geom_point(data = slice_sample(combined, n = 2000), aes(x = x, y = y), col = "grey") +
  geom_hex(data = mcmc, aes(x = x, y = y, fill = ..density..), bins = 100) +
  geom_point(data = sp, aes(x = x, y = y), col = "darkcyan", fill = "cyan2", shape = 21, size = 2) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~ chain, nrow = 5) + 
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())

# X coordinate density
dx <- ggplot(combined, aes(x = x, fill = type, col = type)) + 
  geom_density(alpha = 0.3) +
  theme_bw() +
  scale_y_reverse() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
# Y coordinate density
dy <- ggplot(combined, aes(x = y, fill = type, col = type)) + 
  geom_density(alpha = 0.3) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
# XY overlay
dxy <- ggplot(combined_sub, aes(x = x, y = y, col = type)) +
  geom_point(size = 1, alpha = 0.2) +
  theme_bw()

combo <- dxy + dy + dx + guide_area() + plot_layout(guides = 'collect')
sensor_plot <- dgrid | combo

sensor_plot



ggplot(mcmc, aes(x = x)) + 
  geom_density(alpha = 0.3)

dx






# Plotting function
p <- function(d, n, x, y, lims){
  plot_cor(d, paste0("x", n), paste0("y", n)) +
    geom_vline(xintercept = x, lty = 2) +
    geom_hline(yintercept = y, lty = 2) +
    xlim(lims[1], lims[2]) +
    ylim(lims[3], lims[4]) +
    theme(legend.position = "none")
}

# Get standard plot limits
limits <- sensor_output_tempered$output
limits <- limits[round(seq(1, nrow(limits), length.out = 2000)), ]
limits <- limits[limits$phase == "sampling", 6:13]
limits <- apply(limits, 2, range)
limits[1,] <- limits[1,] - 0.1
limits[2,] <- limits[2,] + 0.1

# Plots
p1 <- p(sensor_output, 1, 0.57, 0.91, c(-0.2, 0.8, 0.3, 1.2))
p2 <- p(sensor_output, 2, 0.1, 0.37, c(-0.25, 1.25, -0.5, 1))
p3 <- p(sensor_output, 3, 0.26, 0.14, c(0, 1, -0.1, 0.75))
p4 <- p(sensor_output, 4, 0.85, 0.04, c(-1.2, 1.75, -0.75, 2))
p5 <- p(sensor_output_tempered, 1, 0.57, 0.91,  c(-0.2, 0.8, 0.3, 1.2))
p6 <- p(sensor_output_tempered, 2, 0.1, 0.37, c(-0.25, 1.25, -0.5, 1))
p7 <- p(sensor_output_tempered, 3, 0.26, 0.14, c(0, 1, -0.1, 0.75))
p8 <- p(sensor_output_tempered, 4, 0.85, 0.04, c(-1.2, 1.75, -0.75, 2))
# Combine plots
sensors <- (p1 / p2 / p3 / p4) | (p5 / p6 / p7 / p8)

# save to file
ggsave("Figure_sensors/sensors.png", sensors, width = 6, height = 8)
