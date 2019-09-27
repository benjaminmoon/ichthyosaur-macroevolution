#!/usr/bin/env Rscript

# Functions for completing rates analyses.
# From Moon & Stubbs <Title>, <doi>

# Some original scripts/functions from Close _et al._ (2015, _Curr. Biol.,_
# doi:10.1016/j.cub.2015.06.047). Supplementary information found at
# http:// doi.org/10.5061/dryad.760sc.

# Contents
#   - divTaxPhy
#   - fillInsignificant
#   - fillMissingBins

library(paleotree)
library(pbapply)

divTaxPhy <- function (bins, ages, time_trees, clus = NULL) {
  # Calculates binned taxic and phylogenetic diversity from first and last
  # occurrences and time-scaled trees. Progress bars are shown for when using
  # large, or many trees; can also be parallelised by givng the cluster name for
  # the environment. Requires `pbapply` and `paleotree` packages.
  #
  # Args:
  #   bins: a list of data frames givng bin maximum and minimum ages, with names
  #         (if desired) as row names.
  #   ages: a data frame with taxon names as row names and giving first and last
  #         occurrences.
  #   time_trees: a multiPhylo object of time-scaled trees.
  #   clus: the name of a cluster object created by makeCluster to run in
  #         parallel upon.
  #
  # Returns:
  #   A list of data frames (one for each binning scheme) with columns for bin
  #   start and end dates, per-bin taxic diversity, per-bin mean and median
  #   phylogenetic diversity and total ranges based on supplied trees.
  div_binlist <- lapply(bins, function (bin_scheme) {
    # phylogenetic diversity using paleotree
    # keeps only the diversity (column 3)
    phydiv <- pblapply(cl = clus, time_trees, function (tree) {
      paleotree::phyloDiv(tree, int.times = bin_scheme, plot = FALSE)[, 3]
    })
    # combine into table
    phydiv <- cbind(bin_scheme,
                    matrix(unlist(phydiv),
                           nrow  = nrow(bin_scheme),
                           byrow = FALSE))
    # calculate mean, median, and range
    data_summ <- list(median, mean, max, min)
    summ_phydiv <- lapply(data_summ, function (summ) {
      apply(phydiv[, 3:ncol(phydiv)], 1, summ, na.rm = FALSE)
    })

    # Create data.frame of values for plotting: taxic diversity, mean and median
    # phylogenetic diversity etc.
    div_df <- data.frame(max_age = phydiv[, 1],
                         min_age = phydiv[, 2],
                         mid_age = rowMeans(phydiv[, 1:2]),
                         tax_div = taxicDivCont(ages,
                                                int.times = bin_scheme,
                                                plot      = FALSE)[, 3],
                         phy_meddiv  = summ_phydiv[[1]],
                         phy_meandiv = summ_phydiv[[2]],
                         phy_maxdiv  = summ_phydiv[[3]],
                         phy_mindiv  = summ_phydiv[[4]])

    # return data.frame
    return(div_df)
  })
}

fillInsignificant <- function (X) {
  # Removes NULL runs and replaces insignificant runs (organising for plotting)
  # following Claddis::DiscreteCharacterRate runs on multiple trees. From Close
  # _et al._ 2015.
  #
  # Args:
  #   X: the list of Claddis::DiscreteCharacterRates runs
  #
  # Returns:
  #   A list of Claddis::DiscreteCharacterRate runs with appropriate
  #   replacements.
  # remove NULL values from the resultant list
  X[sapply(X, is.null)] <- NULL

  # fill insignificant runs with per.bin.rate matrices
  if (!is.null(X[[1]]$per.bin.rates.ib)) {
    X.ib <- X; for (i in 1:length(X.ib)) {
      X.ib[[i]]$per.bin.rates <- X.ib[[i]]$per.bin.rates.ib
      # duplicate object with only internal branch comparisons
      # for ease of plotting
    }
    X.tb <- X; for (i in 1:length(X.tb)) {
      X.tb[[i]]$per.bin.rates <- X.tb[[i]]$per.bin.rates.tb
      # duplicate object with only terminal branch comparisons
      # for ease of plotting
    }
  }

  # fill missing bins
  X <- fillMissingBins(X)
  if (exists("X.ib")) {
    X.ib  <- fillMissingBins(X.ib)
    X.tb  <- fillMissingBins(X.tb)
  }
}


fillMissingBins <- function (X) {
  # "Here follows some disgustingly convoluted code to convert the single
  # per.bin.rate value output by insignificant runs into a matrix akin to those
  # from significant runs." From Close _et al._ 2015.
  sig.runs <- which(sapply(X, FUN = function(x) is.matrix(x$per.bin.rates))) # returns vector of indices identifying significant runs of DiscCharRates in X
  insig.runs <- which(sapply(X, FUN = function(x) !is.matrix(x$per.bin.rates))) # vector of indices identifying insignificant runs
  sig.nrows <- nrow(X[[sig.runs[1]]]$per.bin.rates) # how many rows do significant $per.bin.rates instances contain?

  if (length(insig.runs) > 0) {
    # the following loop fills in $per.branch.rates for insigificant runs with a matrix like those for significant runs, allowing the plotting code to proceed without errors
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$per.bin.rates <- as.numeric(strsplit(strsplit(X[[insig.runs[i]]][['per.bin.rates']], "of ")[[1]][3], "is preferred.")[[1]][1])
    }
    
    # for each insigificant run in X.branchpies copy a significant $per.bin.rates matrix into a new tmp object
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates <- X[[sig.runs[1]]]$per.bin.rates
    }
    
    # populate "in.rate" column with single rate estimate
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates[,"in.rate"] <- rep(X[[insig.runs[i]]]$per.bin.rates, sig.nrows)
    }
    
    # same for "out.rate"
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates[,"out.rate"] <- rep(X[[insig.runs[i]]]$per.bin.rates, sig.nrows)
    }
    
    # fill "ml.signif.hi" with zeroes to reflect non-significant results
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates[,"ml.signif.hi"] <- rep(0, sig.nrows)
    }
    
    # fill "ml.signif.lo" with zeroes
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates[,"ml.signif.lo"] <- rep(0, sig.nrows)
    }
    
    # fill "ml.signif" with zeroes
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates[,"ml.signif"] <- rep(0, sig.nrows)
    }
    
    # copy newly-constructed temporary $per.bin.rates matrices into proper location
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$per.bin.rates <- X[[insig.runs[i]]]$tmp.per.bin.rates
    }
    
    # delete temporary matrix
    for (i in 1:length(insig.runs)) {
      X[[insig.runs[i]]]$tmp.per.bin.rates <- NULL
    }  
  }

  # return data
  return(X)
}