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

# Plot Metropolis coupling acceptance rates
plot_mc_acceptance(sensor_output_tempered)

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
limits <- limits[limits$phase == "sampling" & limits$rung == 10, 7:14]
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
