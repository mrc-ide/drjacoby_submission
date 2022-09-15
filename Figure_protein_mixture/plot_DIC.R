# plot_comp_means
#
# Author: Bob Verity
# Date: 2021-10-05
#
# Purpose:
# Read in MCMC output from run_TAGM.R and produce ggplot of DIC over a range of
# K for each of the MCMC sampling methods.
#
# ------------------------------------------------------------------

library(magrittr)
library(ggplot2)

# ------------------------------------------------------------------

# define plotting parameters
col_vec <- c("#cc282b", "#53b3cb")

# read in MCMC and extract DIC values
K_vec <- 2:15
df_res <- data.frame(K = K_vec,
                     DIC_standard = NA,
                     DIC_PT = NA)
for (i in seq_along(K_vec)) {
  mcmc_basic <- readRDS(sprintf("/Users/rverity/Desktop/protein_mixture_output/TAGM_basic_%s", K_vec[i]))
  mcmc_rungs <- readRDS(sprintf("/Users/rverity/Desktop/protein_mixture_output/TAGM_rungs_%s", K_vec[i]))
  df_res$DIC_standard[i] <- mcmc_basic$diagnostics$DIC_Gelman
  df_res$DIC_PT[i] <- mcmc_rungs$diagnostics$DIC_Gelman
}

# get minimum values by each method
K_standard <- df_res$K[which.min(df_res$DIC_standard)]
K_PT <- df_res$K[which.min(df_res$DIC_PT)]

# plot
plot1 <- df_res %>%
  dplyr::filter(K > 1) %>%
  dplyr::rename(Standard = DIC_standard,
                PT = DIC_PT) %>%
  tidyr::pivot_longer(cols = -K) %>%
  dplyr::mutate(name = c("Basic MH", "PT")[match(name, c("Standard", "PT"))]) %>%
  dplyr::mutate(name = factor(name, levels = c("Basic MH", "PT"))) %>%
  ggplot(aes(x = K, y = value, col = name)) + theme_bw() +
  geom_line() + geom_point() +
  scale_color_manual(values = col_vec) +
  facet_wrap(~name, scales = "free", ncol = 1) +
  scale_x_continuous(breaks = K_vec) +
  xlab("Number of mixture components (K)") + ylab("Deviance Information Criterion") +
  theme(panel.spacing.x = unit(1, "cm"),
        strip.text.x = element_text(family = "Arial", size = 12, hjust = 0),
        strip.background.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank())

plot1

# save plot grob to file
saveRDS(plot1, file = "Figure_protein_mixture/output/protein_mixture_DIC_plot.rds")

# save image to file
ggsave("Figure_protein_mixture/protein_mixture_DIC.png", plot1, width = 4, height = 6)
