#!/usr/bin/env Rscript

# 3. Discrete rates analyses #
# ========================== #

# Script from Moon & Stubbs Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

# required libraries
library(Claddis)
library(doParallel)
library(magrittr)
library(parallel)
library(pbapply)
library(rlist)

# required scripts
source("functions/plotting_functions.R")
source("functions/rates_figures.R")
source("functions/rates_functions.R")
source("functions/sendEmailNotification.R")

# for parallelisation
no_cores <- detectCores() - 1
clus <- makeCluster(no_cores, outfile = "")
registerDoParallel(clus)


# 3.1. Gather data #
# ---------------- #

# load tree, age, and cladistic data
nexus_data <- ReadMorphNexus("data/matrix.nex")
bins <- list(epoch    = read.table("data/epoch_bins.tsv",
                                   sep       = "\t",
                                   header    = TRUE,
                                   row.names = 1),
             tenMa_bJ = data.frame(max_age = seq(251.3, 100, -10),
                                   min_age = seq(241.3, 90, -10)))
taxon_ages <- read.table("./data/ichthyosaur_occurrences.tsv",
                         header    = TRUE,
                         row.names = 1,
                         sep       = "\t")
time_trees <- readRDS("output/time_trees.rds")


# 3.2. Diversity curves #
# --------------------- #

# calculate diversity using div_taxphy function from Hedman-scaled trees
# (included in R/functions/phylo_div.R)
hedman_diversity <- list.find(time_trees,
                              scaling == "Hedman",
                              n = length(time_trees)) %>%
                    lapply("[[", "tree") %>%
                    divTaxPhy(bins, taxon_ages, time_trees = ., cl = clus)
saveRDS(hedman_diversity, file = "output/hedman_diversity.rds")

{
  # plot diversity curves
  cairo_pdf("fig/figS1-hedman_diversity.pdf",
            width  = 5,
            height = 7)
  # set-up plot
  par(mfrow   = c(2, 1),
      mar     = c(0, 0, 0, 0),
      oma     = c(5, 3, 0, 0) + 0.25,
      mgp     = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  fig_diversity(hedman_diversity)
  # end plotting to pdf
  dev.off()
}

# calculate diversity from MBL-scaled trees
mbl_diversity <- list.find(time_trees,
                           scaling == "mbl",
                           n = length(time_trees)) %>%
                 lapply("[[", "tree") %>%
                 divTaxPhy(bins, taxon_ages, time_trees = ., cl = clus)
saveRDS(mbl_diversity, file = "output/mbl_diversity.rds")

{
  # plot diversity curves
  cairo_pdf("fig/figS2-mbl_diversity.pdf",
            width  = 5,
            height = 7)
  # set-up plot
  par(mfrow   = c(2, 1),
      mar     = c(0, 0, 0, 0),
      oma     = c(5, 3, 0, 0) + 0.25,
      mgp     = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  fig_diversity(mbl_diversity)
  # end plotting to pdf
  dev.off()
}


# 3.3. Discrete character rates #
# ----------------------------- #

# export packages and variables to cluster for parallelisation
clusterCall(clus, function () {
  # required libraries
  library(Claddis)
  library(magrittr)
})

# calculate rates on trees (added code forms a progress bar)
# null values are included in the resultant list; to be removed later
discrete_rates <- lapply(seq_along(bins), function (scheme) {
                    # get bin intervals
                    intervals <- c(bins[[scheme]][, 1],
                                   min(bins[[scheme]]))

                    # apply to each `time_trees`
                    lapply(time_trees, function (tree) {
                      # combine into list
                      return(list(nexus   = nexus_data,
                                  tree    = tree$tree,
                                  scaling = tree$scaling,
                                  bins    = list(bins = intervals,
                                                 name = names(bins[scheme]))))
                    })
                  }) %>%
                  unlist(recursive = FALSE) %>%
                  pblapply(cl = clus, function (run) {
                    # run discrete rates analysis, return NULL if failed
                    r <- DiscreteCharacterRate(run$tree,
                                               run$nexus,
                                               time.bins = run$bins$bins,
                                               alpha     = 0.01) %>%
                         tryCatch(error = function (e) NULL)

                    # return rates
                    return(c(r,
                             bins    = run$bins$name,
                             scaling = run$scaling))
                  })
# send notification email when finished
sendEmailNotification(subject = "Notification: DiscreteCharacterRate success",
                      body    = "DiscreteCharacterRate analyses completed successfully. :)")

# remove failed attempts and organise for plotting
# taken from Close _et al._ (2015)
discrete_rates <- fillInsignificant(discrete_rates)
saveRDS(discrete_rates, file = "output/discrete_rates.rds")

# spaghetti plotting: plot time series of mean rates across all trees
{
  cairo_pdf("fig/fig3-rates_spaghetti.pdf",
            width = 8,
            height = 7)
  # set-up plot
  par(mfrow   = c(2, 1),
      mar     = c(0, 0, 0, 0),
      oma     = c(5, 3, 0.25, 0) + 0.25,
      mgp     = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  list.find(discrete_rates,
            scaling == "Hedman",
            n = length(discrete_rates)) %>%
  fig_spaghettiHedman()
  # end plotting to pdf
  dev.off()
}

# spaghetti plot of MBL trees
{
  cairo_pdf("fig/figS8-rates_MBLspaghetti.pdf",
            width = 8,
            height = 7)
  # set-up plot
  par(mfrow   = c(2, 1),
      mar     = c(0, 0, 0, 0),
      oma     = c(5, 3, 0, 0) + 0.25,
      mgp     = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  list.find(discrete_rates,
            scaling == "mbl",
            n = length(discrete_rates)) %>%
  fig_spaghettiHedman()
  # end plotting to pdf
  dev.off()
}

# stop cluster
stopCluster(clus)
