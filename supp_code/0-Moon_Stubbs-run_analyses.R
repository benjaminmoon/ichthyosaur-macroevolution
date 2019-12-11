#!/usr/bin/env Rscript

# The R scripts sourced here include the analyses for Moon & Stubbs 'Early high
# disparity and rates in the evolution of ichthyosaurs', Communications
# Biology

# Installation of required packages #
# --------------------------------- #

# packages on CRAN
list_of_packages <- c("ape", "colorspace", "devtools", "dispRity", "doParallel",
                      "foreach", "ggplot2", "keyringr", "magrittr", "mailR",
                      "paleotree", "paleoTS", "pbapply", "rlist", "strap",
                      "vegan")
new_packages <- list_of_packages[!(list_of_packages %in%
                                   installed.packages()[, "Package"])]

# change CRAN repository as desired
if (length(new_packages)) {
  install.packages(new_packages, repos = "https://www.stats.bris.ac.uk/R/")
}

# recent package versions from GitHub
# using phytools package from Dec 2017 to make version available with Claddis
# below
if (!require(phytools)) devtools::install_github("liamrevell/phytools")
# this code was created before recent changes to the Claddis package (v0.3),
# therefore we specify an older version to install and run
if (packageVersion("Claddis") != "0.2") {
  # Downgrade Claddis to suitable version
  remove.packages("Claddis")
}
if (!require(Claddis)) devtools::install_github("graemetlloyd/Claddis",
                                                ref = "4f7f1bf")

# check for directories
directories <- c("./fig", "./output", "./tbl")
missing_dir <- directories[!(directories %in% list.dirs())]

# create new directories
for (dir in missing_dir) dir.create(dir)

# load script
source("functions/sendEmailNotification.R")


# 1. Disparity analyses #
# ===================== #

# Conversion of cladistic data to taxon-distance matrices. Pairwise distances
# on distance matrices, sum of ranges, sum of variances, and centroid distance
# on PCo analyses of the distance matrices, with and without Caillez (1983)
# negative eigenvalue correction. All on epoch-length and 10 Ma bins.

source("1A-disparity_analyses.R")

# Statistical comparisons of the disparity data: PERMANOVA on morphospace and
# t-tests on disparity.

source("1B-disparity_stats.R")


# 2. Time-scaling trees #
# ===================== #

# Time-scaling of trees from Moon (2017, _J. Syst. Palaeontol.,_
# doi:10.1080/14772019.2017.1394922) using Hedman (2010) and minimum branch
# length methods.

source("2-time_scaling.R")


# 3. Discrete rates analyses #
# ============================

# Discrete character rates analyses in epoch-length and 10 Ma bins.

source("3-rates_analyses.r")


# send notification email when finished
sendEmailNotification(subject = "Ichthyosaur analyses completed",
                      body    = "Woohoo! All the analyses have finished.\nLet's rock to the pub.")
