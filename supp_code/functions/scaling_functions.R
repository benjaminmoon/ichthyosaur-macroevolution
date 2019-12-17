#!/usr/bin/env Rscript

# Functions for time-scaling trees.
# From Moon & Stubbs "Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

# Contents
#   - timeHedman

timeHedman <- function (ages) {
  # Samples taxon occurrences from a uniform distribution then formats them for
  # use in Hedman (2010) tree time-scaling.
  #
  # Args:
  #   ages: two-column matrix of per-taxon first (FAD) and last (LAD) occurrence
  #         dates with taxon row names
  #
  # Returns:
  #   A vector of sampled per-taxon occurrence dates with row names.
  # sample occurrence dates
  d_samp <- ages %>%
            apply(1, function (taxon) runif(1, min = taxon[2], max = taxon[1]))

  # return named vector
  return(d_samp)
}
