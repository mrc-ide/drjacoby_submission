
# process_data.R
#
# Author: Bob Verity
# Date: 2021-01-14
#
# Purpose:
# Load raw data and get into format for subsequent analysis.
#
# ----------------------------------------------------------------
# install and load packages

#install.packages("BiocManager")
#BiocManager::install("pRolocdata")
library(pRolocdata)

#BiocManager::install("pRoloc")
library(pRoloc)

# ----------------------------------------------------------------

# load raw data from pRolocdata package
#data("hyperLOPIT2015")
dim(hyperLOPIT2015)

# convert to matrix
# method taken from https://github.com/lgatto/pRoloc/blob/02486222d9874c29d0fd8fb9fa207c7993546280/R/plotting-ellipse.R
markersubset <- pRoloc::markerMSnSet(hyperLOPIT2015)
mydata <- Biobase::exprs(markersubset)

# save to file. This dataset is large so processed data will be kept out of the
# Git repos. Other users will need to alter this filepath manually
saveRDS(mydata, file = "/Users/rverity/Desktop/hyperLOPIT2015_processed.rds")
