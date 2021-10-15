# plot_comp_means
#
# Author: Bob Verity
# Date: 2021-10-05
#
# Purpose:
# Read in MCMC output from run_TAGM.R and produce ggplot of estimated
# component means vs. data.
#
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ------------------------------------------------------------------

# read in drjacoby project
K <- 8
mcmc_basic <- readRDS(sprintf("/Users/rverity/Desktop/protein_mixture_output/TAGM_basic_%s", K))
mcmc_rungs <- readRDS(sprintf("/Users/rverity/Desktop/protein_mixture_output/TAGM_rungs_%s", K))

# extract data
dat <- mcmc_basic$parameters$data

# combine posterior draws from both MCMC runs
mcmc_basic$output$model = "basic"
mcmc_rungs$output$model = "rungs"
output <- rbind(mcmc_basic$output,
                mcmc_rungs$output)

# subset to sampling iterations
mcmc_samples <- output %>%
  dplyr::filter(phase == "sampling") %>%
  dplyr::filter(chain <= 3)

# get parameters into long format
param_names <- c(sprintf("mu%s_x", 1:K), sprintf("mu%s_y", 1:K))
mcmc_long <- mcmc_samples %>%
  dplyr::select(model, chain, iteration, param_names) %>%
  tidyr::pivot_longer(param_names)

# group all mu into x vs y coordinates
mcmc_long$name <- substr(mcmc_long$name, nchar(mcmc_long$name), nchar(mcmc_long$name))

# split x and y back into separate columns
tmp <- mcmc_long %>%
  dplyr::filter(name == "x")
df_plot <- data.frame(model = tmp$model,
                      chain = tmp$chain,
                      x = tmp$value,
                      y = subset(mcmc_long, name == "y")$value)

# format plotting names
df_plot$model_name <- c("Standard", "PT")[match(df_plot$model, c("basic", "rungs"))]
df_plot$model_name <- factor(df_plot$model_name, levels = c("Standard", "PT"))
df_plot$chain_name <- sprintf("Chain %s", df_plot$chain)

# plot
plot1 <- ggplot() + theme_bw() +
  geom_point(aes(x = x, y = y, col = "Data"), size = 0.1, data = as.data.frame(dat)) +
  geom_point(aes(x = x, y = y, col = "Posterior draws of\ncomponent means"), size = 0.1, data = df_plot) +
  scale_color_manual(values = c(grey(0.8), "darkblue"), name = "") +
  xlab("PC1 (48% variance)") + ylab("PC2 (23% variance)") +
  facet_grid(rows = vars(chain_name),
             cols = vars(model_name),
             switch = "y") +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        panel.spacing.x = unit(1, "lines"))
plot1

# save plot grob to file
saveRDS(plot1, file = "Figure_protein_mixture/output/protein_mixture_means_plot.rds")

# save image to file
ggsave("Figure_protein_mixture/protein_mixture_means.png", plot1, width = 7, height = 7)
