#!/usr/bin/env Rscript

# 1B. Disparity comparisons #
# ========================= #

# Script from Moon & Stubbs Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

# required libraries
library(doParallel)
library(magrittr)
library(parallel)
library(pbapply)
library(vegan)

# required scripts
source("functions/disparity_stats_functions.R")

# for parallel processing
# register parallel cores for disparity loop below
# change number of cores as available/desired
no_cores <- detectCores() - 1
clus <- makeCluster(no_cores, outfile = "")
registerDoParallel(clus)

# get required data
bin_data <- readRDS("output/bin_data.rds")
disparity <- readRDS("output/disparity.rds")
pco_data <- readRDS("output/pco_data.rds")

# 1.6. Pairwise morphospace PERMANOVA #
# ----------------------------------- #

# add epoch bins to pco_data data.frames
permanova_data <- lapply(pco_data, function (pco) {
  # sort pco data into bins (minus outgroup)
  pc_dat <- list()
  for (bin in seq_along(bin_data$epochs$bin_data[-9])) {
    # get taxa in the bin
    bin_tax <- bin_data$epochs$bin_data[[bin]]

    # get pc data for taxa in the bin and add group name
    pc_dat[[bin]] <- as.data.frame(pco$pco[which(rownames(pco$pco) %in% bin_tax), ])

    # add group column
    pc_dat[[bin]]$bin_nam <- names(bin_data$epochs$bin_data[bin])
  }
  pc_dat <- do.call("rbind", pc_dat)

  # return list
  return(list(dist = pco$dist,
              corr = pco$corr,
              pco  = pc_dat))
})

# export function
clusterExport(clus, c("pairwise_adonis", "adonis"))

# pairwise PERMANOVA contrasts
permanova <- pblapply(permanova_data, function (dat) {
  # number of columns to use (not bin column)
  cols <- 1:(ncol(dat$pco) - 1)

  # permanova test
  out <- pairwise_adonis(dat$pco[, cols],
                         factors = dat$pco[, "bin_nam"],
                         sim_method = "euclidean",
                         p_adjust_m = "fdr")

  # return list
  return(list(dist      = dat$dist,
              corr      = dat$corr,
              permanova = out))
}, cl = clus)

# export as CSV files
for (tab in permanova) {
  # data
  dat <- tab$permanova

  # name
  nam <- paste0("tbl/permanova-", tab$dist, "-", tab$corr, ".csv")

  # write TSV
  write.csv(dat, file = nam)
}


# 1.7. Pairwise disparity t-tests #
# ------------------------------- #

# perform pairwise t-tests on disparity between each bin
ttest_res <- pblapply(disparity, function (run) {
               # pairwise t-test with false discovery rate adjustment
               ttest_dat <- pairwise_ttest(run$summ_test, p_adjust_m = "fdr")

               # return data
               return(list(dist = run$dist,
                           corr = run$corr,
                           bins = run$bins,
                           metric = run$metric,
                           ttest = ttest_dat))
})

# export t-tests as CSV files
for (tab in ttest_res) {
  # data
  dat <- tab$ttest

  # name
  nam <- paste0("tbl/ttest-",
                tab$dist,
                "-",
                tab$corr,
                "-",
                tab$bins,
                "-",
                tab$metric,
                ".csv")

  # write CSV
  write.csv(dat, file = nam, row.names = TRUE)
}

# stop cluster
stopCluster(clus)
