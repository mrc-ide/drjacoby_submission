# plot_combined
#
# Author: Bob Verity
# Date: 2021-10-04
#
# Purpose:
# Read in protein mixture plotting graphical objects and produce combined plot.
#
# ------------------------------------------------------------------

# read in plotting grobs
plot1 <- readRDS("Figure_protein_mixture/output/protein_mixture_means_plot.rds") +
  labs(tag = "A")
plot2 <- readRDS("Figure_protein_mixture/output/protein_mixture_DIC_plot.rds") +
  labs(tag = "B")

# produce combined plot
b <- 6
m <- cbind(matrix(1, nrow = b, ncol = 4), matrix(c(3, rep(2, b - 2), 3), nrow = b, ncol = 3))
plot_c <- gridExtra::grid.arrange(plot1, plot2, layout_matrix = m)

# save image to file
ggsave("Figure_protein_mixture/protein_mixture_combined.png", plot_c, width = 8, height = 6)
