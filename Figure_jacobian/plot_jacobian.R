# plot_jacobean.R
#
# Author: Pete Winskill
# Date: 2021-09-07
#
# Purpose:
# Produce schematic plot demonstrating the logic behind the Jacobean
# transformation.
#
# ------------------------------------------------------------------

### Load packages ##############################################################
library(ggplot2)
library(dplyr)
library(patchwork)
################################################################################

### Functions for transformations ##############################################
get_phi <- function(x, a = 0){
  log(x - a)
}
get_theta <- function(x, a = 0){
  exp(x) + a
}
################################################################################

### Plot 1: Original range #####################################################
d1 <- data.frame(x = 10, y = 0, xmin = 0, xmax = 30)

p1 <- ggplot(data = d1, aes(x = x, y = y, xmin = xmin, xmax = xmax)) +
  geom_point(size = 3, col = "darkblue") +
  geom_point(x = 0, y = 0, size = 3, shape = "|", col = "darkblue") +
  geom_linerange(size = 1, col = "darkblue") + 
  geom_text(x = 30, y = -0.02, label = expression(infinity)) +
  geom_text(x = 0, y = -0.02, label = "0") +
  geom_segment(aes(x = 25, xend = 30, y = 0, yend = 0),  arrow = arrow(length = unit(0.2, "cm")), col = "darkblue", size = 1) + 
  scale_x_continuous(limits = c(-30, 30)) +
  theme_void() +
  ggtitle("1. Original range")
################################################################################

### Plot 2: Transform and propose ##############################################
# Line and parameter value 
d1_transform <- d1 %>%
  mutate(x = get_phi(x),
         xmin = -30,
         xmax = 30)
# New proposed value
o1 <- data.frame(x = 4.5, y = 0)
# Proposal distribution (normal)
l3 <- data.frame(x = seq(-30, 30, 0.1)) %>%
  mutate(y = dnorm(x, d1_transform$x, 2))
# Arrow from current parameter to proposed
l4 <- data.frame(y = 0 + 0.03, x = o1$x,  xend = d1_transform$x)

p2 <- ggplot() +
  geom_line(dat = l3, aes(x = x, y = y), col = "springgreen4", size = 1) +
  geom_point(aes(x = 0, y = 0), size = 3, shape = "|", col = "darkblue") +
  geom_point(data = d1_transform, aes(x = x, y = y), size = 3, col = "darkblue") +
  geom_linerange(data = d1_transform, aes(x = x, y = y, xmin = xmin, xmax = xmax), size = 1, col = "darkblue") + 
  geom_point(data = o1, aes(x = x, y = y), shape = 18, size = 4, col = "darkred") +
  geom_curve(data = l4, aes(y = y, yend = y, x = x, xend = xend),
             arrow = arrow(length = unit(0.05, "npc"), ends = "first", type = "open"),
             col = "darkred", curvature = 1, size = 0.8) +
  geom_segment(aes(x = 25, xend = 30, y = 0, yend = 0),  arrow = arrow(length = unit(0.2, "cm")), col = "darkblue", size = 1) + 
  geom_segment(aes(x = -25, xend = -30, y = 0, yend = 0),  arrow = arrow(length = unit(0.2, "cm")), col = "darkblue", size = 1) + 
  geom_text(aes(x = 30, y = -0.03), label = expression(infinity)) +
  geom_text(aes(x = -30, y = -0.03), label = expression(-infinity)) +
  geom_text(aes(x = 0, y = -0.02), label = "0") +
  scale_x_continuous(limits = c(-30, 30)) +
  theme_void() +
  ggtitle("2. Transform and propose")
################################################################################

### 3. Backtransform, accept/reject ############################################
#  Back transform the proposed parameter (tweaked for aesthetics)
d1_transform_back <- d1_transform %>%
  mutate(x = 18,
         xmin = 0,
         xmax = 30)
# Back transform the proposal distribution
l3_transform_back <- l3 %>%
  mutate(x = get_theta(x)) %>%
  filter(x <=30, x >= -30)

p3 <- ggplot(data = d1_transform_back, aes(x = x, y = y, xmin = xmin, xmax = xmax)) +
  #geom_line(dat = l3_transform_back, aes(x = x, y = y), col = "springgreen1") +
  geom_point(x = 0, y = 0, size = 3, shape = "|", col = "darkblue") +
  geom_linerange(size = 1, col = "darkblue") + 
  geom_point(shape = 18, size = 4, col = "darkred") +
  geom_segment(aes(x = 25, xend = 30, y = 0, yend = 0),  arrow = arrow(length = unit(0.2, "cm")), col = "darkblue", size = 1) + 
  geom_text(x = 30, y = -0.02, label = expression(infinity)) +
  geom_text(x = 0, y = -0.02, label = "0") +
  scale_x_continuous(limits = c(-30, 30)) +
  theme_void() +
  ggtitle("3. Backtransform, accept/reject")
################################################################################

### Combine and save ###########################################################
p <- (p1 / p2 / p3) + plot_layout(heights = c(1, 1.5, 1))
ggsave("Figure_jacobian/jacobian.png", p, height = 4, width = 7, dpi = 500)
################################################################################
