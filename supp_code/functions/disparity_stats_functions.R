#!/usr/bin/env Rscript

# Functions for script `disparity_stats.R` for performing statistical test on
# disparity analyses.
# From Moon & Stubbs "Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

library(vegan)

# Contents:
#   - pairwise_adonis
#   - pairwise_ttest
#   - perBinReplicates
#   - t_test

pairwise_adonis <- function (x, factors, sim_method = "euclidean",
                             p_adjust_m = "fdr", perm = 9999) {
  # script created by Pedro Arbizu; downloaded from
  # https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
  #
  # Calculates PERMANOVA contrasts on a set of PC(o) data matrix with a defined
  # grouping (e.g. binning).
  #
  # Args:
  #   x: the data matrix of taxa and observations to use for contrasts.
  #   factors: the grouping to use, can be within the matrix, or separately, as
  #            long as it is matched to the matrix.
  #   sim_method: simulation method for data
  #   p_adjust_m: p-value adjustment method
  #
  # Returns:
  #   Data frame of pairwise PERMANOVA F-statistics, p-values, and adjusted p-values.
  #   pairwise comparisons.
  co <- combn(unique(factors), 2)

  # set up variables used
  pairs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()

  # pairwise PERMANOVA
  for (elem in 1:ncol(co)){
    # comparable bins
    comp_bin <- c(co[1, elem], co[2, elem])

    # sort the factors/taxa to use
    fact <- factors %in% comp_bin

    # do PERMANOVA
    ad <- vegan::adonis(x[fact, ] ~ factors[fact],
                        method = sim_method, permutations = perm)

    # organise output data
    pairs <- c(pairs, paste(co[1, elem], "vs", co[2, elem]))
    F.Model <- c(F.Model, ad$aov.tab[1, 4])
    R2 <- c(R2, ad$aov.tab[1, 5])
    p.value <- c(p.value, ad$aov.tab[1, 6])
  }

  # adjust p-values
  p.adjusted <- p.adjust(p.value, method = p_adjust_m)

  # output data frame
  pairw.res <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)

  # return data
  return(pairw.res)
}

pairwise_ttest <- function (dat, p_adjust_m = "fdr") {
  # Pairwise t-tests between per-bin disparity from a distance matrix.
  #
  # Args:
  #   dat:        the data matrix to use (include bins here?) with columns
  #               entitled "bs.mean", "st_d", and "samp".
  #   binning:    set of bins to use
  #   replicates: the number of bootstrap replicates to do
  #   method:     the method of disparity estimation
  #   p_adjust_m: p-value adjustment method ("fdr" by default)
  #
  # Returns:
  #   A data.frame of pairwise t-test statistics, and p-values, and adjusted
  #   p-values between bins.
  # combination of groups
  co <- combn(rownames(dat), 2)

  # storage variable
  pairs <- c()
  tt <- list()

  # pairwise t-test
  for (elem in 1:ncol(co)){
    # comparable bins
    comp_bin <- c(co[1, elem], co[2, elem])

    # return NAs on empty bins (too-low sample size)
    na_test <- is.na(c(dat[comp_bin[1], "bs.mean"],
                       dat[comp_bin[2], "bs.mean"]))
    if (any(na_test) == TRUE) {
      tt[[elem]] <- data.frame(t_statistic     = NA,
                               degrees_freedom = NA,
                               p_value         = NA)
    } else {
      # do t-test
      tt[[elem]] <- t_test(dat[comp_bin[1], ], dat[comp_bin[2], ])
    }

    # organise output data
    pairs <- c(pairs, paste(co[1, elem], "vs", co[2, elem]))
  }

  # output data frame
  pairw_res <- do.call(rbind, tt)

  # adjust p-values
  pairw_res$p_adjusted <- p.adjust(pairw_res$p_value,
                                   method = p_adjust_m)

  # replace row names
  rownames(pairw_res) <- pairs

  # return data
  return(pairw_res)
}

perbinReplicates <- function (pco, binning, metric = sum_var,
                              replicates = 500, ci = 95, cluster = clus) {
  # Per-bin bootstrap replicates for several measures of disparity on PC
  # matrices or data.frames. A single PC matrix/data.frame and list of taxon-bin
  # occurrences must be provided. The function subsets taxa into their
  # respective bin, performs the desired disparity calculation for the desired
  # replicates, then returns summary statistics (mean, median, standard
  # deviation) and confidence intervals in a data.frame that can be used for
  # statistical tests. For speed, the replicates have been parallelised.
  #
  # Args:
  #   pco:        the PC matrix or data.frame with taxa (samples) as rows, taxon
  #               names as row names, and PC axis (observations) as columns.
  #   binning:    a list of bins with taxon occurrences in each list item; names
  #               must match the row names of pco.
  #   metric:     one of "sum_rng" (sum of ranges), "sum_var" (sum of
  #               variances), or the same of another function that can be called
  #               upon a matrix or data.frame using `do.call`.
  #   replicates: the number of replicates to perform; default is 10,000.
  #   ci:         confidence interval in integer per cent; default is 95%.
  #   cluster:    the name of the cluster to run on; default is clus.
  #
  # Returns:
  #   A data.frame with rows named as bins and summary data from the replicates:
  #   mean, median, standard deviation, sample size, and upper and lower
  #   confidence interval bounds.
  # divide pc data into bins
  list_pc <- lapply(binning, function (bin) {
               # sub set pc data based on bin-occurrence
               subset(pco, rownames(pco) %in% bin)
             })

  # calculate the confidence interval
  bounds <- (1 - ci / 100) / 2
  bounds <- c(upper = 1 - bounds,
              lower = bounds)

  # perform replicates over list of bin-occurrences
  list_reps <- lapply(list_pc, function (bin_pc) {
                 if (nrow(bin_pc) > 1) {
                   # do replicates; pbsapply give a progress bar for each bin
                   reps <- pbsapply(seq_len(replicates), function (x) {
                             # resample the taxa in the pc data
                             sample <- bin_pc[sample(nrow(bin_pc), replace = TRUE), ]

                             # calculate the metric
                             samp_metr <- metric(sample)

                             # return the data
                             return(samp_metr)
                           }, cl = cluster)

                   # return replicates and sample size
                   return(list(reps = reps,
                               size = nrow(bin_pc)))
                 } else {
                   return(list(reps = NA,
                               size = NA))
                 }
               })

  # summarise list_reps
  reps_summ <- lapply(list_reps, function (bin_reps) {
                 data.frame(mean = mean(bin_reps$reps),
                            median = median(bin_reps$reps),
                            st_d = sd(bin_reps$reps),
                            samp = bin_reps$size,
                            CI_upper = quantile(bin_reps$reps,
                                                probs = bounds[["upper"]],
                                                na.rm = TRUE),
                            CI_lower = quantile(bin_reps$reps,
                                                probs = bounds[["lower"]],
                                                na.rm = TRUE))
               }) %>%
               do.call(rbind, .)
  rownames(reps_summ) <- names(binning)

  # return list of summary data and replicates
  return(list(summary    = reps_summ,
              replicates = list_reps))
}

t_test <- function (x, y) {
  # Calculates Welch's t statistic, which can be used to perform a t-test.
  #
  # Args:
  #   x, y: matrices or data.frames that summarise the sample data. Must include
  #         columns named "bs.mean", "st_d" (standard deviation), and "samp"
  #         (sample size).
  #
  # Returns:
  #   The t-statistic, degrees of freedom, and p-value for the comparison.
  # get difference between the means
  diff_means <- diff(c(x$bs.mean, y$bs.mean))

  # get the variance between the two samples
  varx <- prod(x$samp - 1,
               x$samp,
               x$st_d ^ 2)
  vary <- prod(y$samp - 1,
               y$samp,
               y$st_d ^ 2)
  samp <- sum(x$samp,
              y$samp,
              2)
  var_dist <- (varx + vary) / samp

  # get the sample size differences
  size_diff <- sum(x$samp, y$samp) / prod(x$samp, y$samp)

  # calculate t statistic
  t_stat <- diff_means / sqrt(prod(var_dist, size_diff))

  # calculate degrees of freedom
  deg_freedom <- sum(x$samp,
                     y$samp,
                     -2)

  # calculate significance value
  p_value <- 1 - pt(t_stat, df = deg_freedom)

  # make two-tailed test
  if (p_value > 0.5) {
    p_value <- prod(1 - p_value, 2)
  } else if (p_value < 0.5) {
    p_value <- prod(p_value, 2)
  } else if (p_value == 0.5) {
    p_value <- 1
  }

  # return data
  return(data.frame(t_statistic     = t_stat,
                    degrees_freedom = deg_freedom,
                    p_value         = p_value))
}
