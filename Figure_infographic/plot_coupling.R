
# plot_coupling.R
#
# Author: Bob Verity
# Date: 2020-12-09
#
# Purpose:
# Takes the output of sim_coupling.R and produces a ggplot schematic of how the
# chain moves through temperature rungs.
#
# ------------------------------------------------------------------

# load bobfunctions2 package. If not already installed this can be obtained from
# Github via the command devtools::install_github('bobverity/bobfunctions2')
library(bobfunctions2)

# load other packages
library(ggplot2)

# define double well plotting function
double_well <- function(x, gamma) {
  ret <- exp(-gamma*(x^2 - 1)^2)
  ret / sum(ret)
}

# read in simulation results
sim <- readRDS("Figure_infographic/output/sim_results.rds")
gamma <- sim$gamma
rungs <- sim$rungs
track_full <- sim$track_full

# define plotting parameters
L <- 3
z_max <- 0.07
theta <- 41
phi <- 10
d <- 3

# get curves for each temperature rung
x <- seq(-L, L, l = 1001)
gamma_vec <- seq(gamma, 0, l = rungs)
z <- mapply(function(g) {
  double_well(x, g)
}, gamma_vec)
z[1,] <- z[nrow(z),] <- 0

# get projection matrix
proj_mat <- get_projection(x_lim = c(-L, L),
                           y_lim = c(1, rungs),
                           z_lim = c(0, z_max),
                           theta = theta,
                           phi = phi,
                           d = d)

# project lines
l <- list()
for (i in 1:rungs) {
  df <- project_2d(x, i, z[,i], proj_mat)
  df$rung <- i
  l[[i]] <- df
}
df_curves <- do.call(rbind, l)

# project simulation results
df_track <- project_2d(track_full[,1], track_full[,2], 0, proj_mat)

# produce basic plot
plot1 <- ggplot() + theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())

# add gridlines
grid_x <- seq(-3, 3, 1)
df_gridlines <- expand.grid(c(1,rungs), grid_x)
names(df_gridlines) <- c("y", "x")

df_gridlines_proj <- project_2d(df_gridlines$x, df_gridlines$y, 0, proj_mat)
df_gridlines_proj$group <- rep(seq_along(grid_x), each = 2)

plot1 <- plot1 +
  geom_path(aes(x = x, y = y, group = group), size = 0.1, col = grey(0.6), data = df_gridlines_proj)
plot1


# add x-axis labels
df_xaxis <- data.frame(x = -3:3, y = 0.5)
df_xaxis_proj <- project_2d(df_xaxis$x, df_xaxis$y, 0, proj_mat)
df_xaxis_proj$label = df_xaxis$x

df_xaxis_lab <- project_2d(0.1, 0, 0, proj_mat)

plot1 <- plot1 +
  geom_text(aes(x = x, y = y, label = label), col = grey(0.6), data = df_xaxis_proj) +
  geom_text(aes(x = x, y = y), label = "x", col = grey(0.6), data = df_xaxis_lab)
plot1


# add temperature arrow
delta <- 0.1
df_col <- data.frame(y = seq(1, rungs, delta))
df_col$col <- rev(bobfunctions2::col_hotcold(nrow(df_col)))

df_arrow <- data.frame(x = 3.5, y1 = seq(1, rungs, delta), y2 = seq(1+delta, rungs+delta, delta))
df_arrow_proj <- project_2d(df_arrow$x, df_arrow$y1, 0, proj_mat)
df_arrow_proj2 <- project_2d(df_arrow$x, df_arrow$y2, 0, proj_mat)
df_arrow_proj$x2 <- df_arrow_proj2$x
df_arrow_proj$y2 <- df_arrow_proj2$y
df_arrow_proj$group <- factor(match(df_arrow$y1, df_col$y), levels = 1:nrow(df_col))

plot1 <- plot1 +
  geom_segment(aes(x = x, xend = x2, y = y, yend = y2, col = group), size = 1.5, data = df_arrow_proj) +
  geom_segment(aes(x = x2, xend = x, y = y2, yend = y, col = group), size = 1.5,
               arrow = arrow(length = unit(0.25, "cm")), data = df_arrow_proj[1,]) +
  geom_segment(aes(x = x, xend = x2, y = y, yend = y2, col = group), size = 1.5,
               arrow = arrow(length = unit(0.25, "cm")), data = df_arrow_proj[nrow(df_arrow_proj),]) +
  scale_color_manual(values = df_col$col, guide = "none")
plot1


# add "hot" and "cold" text
ang <- 0
df_hotcold_lab <- data.frame(x = c(4, 4), y = c(1, rungs))
df_hotcold_lab_proj <- project_2d(df_hotcold_lab$x,
                                  df_hotcold_lab$y,
                                  0, proj_mat)
df_hotcold_lab_proj$label <- c("Cold", "Hot")
df_hotcold_lab_proj$group <- factor(match(df_hotcold_lab$y, df_col$y))

plot1 <- plot1 +
  geom_text(aes(x = x, y = y, label = label, col = group), angle = ang, data = df_hotcold_lab_proj)
plot1


# add temperature text
df_yaxis_proj <- project_2d(4.5, mean(1:rungs), 0, proj_mat)

plot1 <- plot1 +
  geom_text(aes(x = x, y = y, label = "Temperature"), angle = ang, size = 5, data = df_yaxis_proj)
plot1


# add curves
df_curves$group <- factor(match(df_curves$rung, df_col$y))

plot1 <- plot1 +
  geom_path(aes(x = x, y = y), data = df_track) +
  geom_point(aes(x = x, y = y), size = 1.5, data = df_track) +
  geom_line(aes(x = x, y = y, col = group, group = group), data = df_curves) +
  geom_polygon(aes(x = x, y = y, fill = as.factor(rung)),
               col = NA, alpha = 0.3, data = df_curves) +
  scale_fill_manual(values = df_col$col[match(1:rungs, df_col$y)], guide = "none")
plot1

# save to file
ggsave("Figure_infographic/infographic_coupling.png")
