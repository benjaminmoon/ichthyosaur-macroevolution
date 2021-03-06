#!/usr/bin/env Rscript

# Scripts for computing disparity.
# From Moon & Stubbs "Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

# Original pairwise distances functions from Close _et al._ (2015, _Curr. Biol.,_
# doi:10.1016/j.cub.2015.06.047). Supplementary information found at
# http://doi.org/10.5061/dryad.760sc.

# Calculate mean and weighted mean pairwise distances from a symmetrical
# dissimilarity matrix, with comparable characters, as produced from
# `Claddis::MorphDistMatrixFast`.

library(dispRity)
library(doParallel)
library(magrittr)

# Contents:
#   - upperTriangle
#   - bootstrapMPD
#   - bootstrapWMPD
#   - doPD
#   - pcCorrelation
#   - rarefyPD
#   - scree
#   - subsetBootstrap

bootstrapMPD <- function(dissim) {
  # Computes the mean pairwise distances between values in a symmetrical dissimilarity matrix.
  # Also bootstraps 10 000 times to get 95% confidence intervals.
  #
  # Args:
  #   dissim: a symmetrical dissimilarity matrix
  #
  # Returns:
  #   A three-column matrix of mean pairwise distance with upper and lower 95% confidence intervals.
  dissim <- upperTriangle(dissim)
  mean <- mean(dissim, na.rm = TRUE)
  Z <- length(dissim[complete.cases(dissim)])
  boot_mean <- vector()
  for (i in 1:500) {
    boot_mean[i] <- mean(dissim[complete.cases(dissim)][sample.int(Z, Z, replace = TRUE)])
  }
  #Lower 0.025 for the mean
  lower <- sort(boot_mean)[length(boot_mean) * 0.025]
  #Upper 0.975 for the mean
  upper <- sort(boot_mean)[length(boot_mean) * 0.975]
  return(cbind(mean, lower, upper))
}

bootstrapWMPD <- function(dissim, compars) {
  # Computes the mean pairwise distances between values in a symmetrical dissimilarity matrix weighted
  # by the proportion of comparable characters. Also bootstraps 500 times to get 95% confidence intervals.
  #
  # Args:
  #   dissim: a symmetrical dissimilarity matrix
  #   compars: a symmetrical matrix of comparable characters
  #
  # Returns:
  #   A three-column matrix of weighted mean pairwise distance with upper and lower 95% confidence intervals.
  dissim <- upperTriangle(dissim)
  compars <- upperTriangle(compars)
  dat <- as.data.frame(cbind(dissim, compars))
  weighted_mean <- sum(dissim * compars, na.rm = TRUE) / sum(compars, na.rm = TRUE)
  Z <- length(dissim[complete.cases(dissim)])
  boot_mean <- vector()
  for (i in 1:500) {
    temp_dat <- dat[complete.cases(dissim), ][sample.int(Z, Z,
      replace = TRUE), ]
    boot_mean[i] <- sum(temp_dat$dissim * temp_dat$compars,
      na.rm = TRUE) / sum(temp_dat$compars, na.rm = TRUE)
  }
  #Lower 0.025 for the mean
  lower <- sort(boot_mean)[length(boot_mean) * 0.025]
  #Upper 0.975 for the mean
  upper <- sort(boot_mean)[length(boot_mean) * 0.975]
  return(cbind(weighted_mean, lower, upper))
}

doPD <- function (binning, dissim, compars) {
  # Performs MPD and WMPD functions on binned occurrence data and outputs a data frame.
  #
  # Args:
  #   binning: a list of bins with taxon assignments matching row and column
  #            names in `dissim` and `compars`
  #   dissim: a symmetrical dissimilarity matrix
  #   compars: a symmetrical matrix of comparable characters
  #
  # Returns:
  #   A list of two elements each with data frames of mean pairwise distances
  #   and weighted mean pairwise distances
  #   respectively for each bin, and confidence intervals.
  # pairwise distances
  mpd <- list()
  mpd <- foreach (i = seq_along(binning)) %dopar% {
    bootstrapMPD(dissim[binning[[i]], binning[[i]]])
  }
  mpd <- do.call(rbind, mpd) %>%
         as.data.frame(stringsAsFactors = FALSE)
  rownames(mpd) <- names(binning)
  # weighted pairwise distances
  wmpd <- list()
  wmpd <- foreach (i = seq_along(binning)) %dopar% {
    bootstrapWMPD(dissim[binning[[i]], binning[[i]]],
                  compars[binning[[i]], binning[[i]]])
  }
  wmpd <- do.call(rbind, wmpd) %>%
          as.data.frame(stringsAsFactors = FALSE)
  rownames(wmpd) <- names(binning)
  # return values
  list(MPD = mpd, WMPD = wmpd) %>% return
}

pcCorrelation <- function (pco, dist) {
  # Correlation tests between a symmetrical distance matrix and the resulting
  # PCo axes by increasing number of PCo axes.
  #
  # Args:
  #   pco:  matrix of taxa (rows) by PCo co-ordinates (columns).
  #   dist: symmetrical distance matrix corresponding to pco.
  #
  # Returns:
  #   A vector of correlation values for each of increasing number of PCo axes.
  # squared distances & PCo
  dist_sq <- dist[lower.tri(dist)] ^ 2

  # correlate to increasing number of PCo axes
  pc_cor <- sapply(seq_along(1:ncol(pco)), function (axes) {
              # PCo squared distances
              pc_dist <- dist(pco[, 1:axes]) ^ 2

              # correlation test
              pc_cor <- cor(dist_sq, pc_dist)

              # return correlation
              return(pc_cor)
            })
}

rarefyPD <- function (d_dist) {
  # Rarefied pairwise distances calculations
  #
  # Args:
  #   d_dist: a list of distance data with associated list of binning data under
  #           ddist and dbin item respectively
  #
  # Returns:
  #   A list of rarefied data frames for each bin and each distacne matrix.
  # Rarefy data
  pblapply(d_dist, cl = clus, function (run) {
    brare <- lapply(run$dbin, function (bin) {
      # get sequence for rarefaction values
      if (length(bin) > 1) {
        nrar <- seq(2, length(bin))
      } else {
        nrar <- length (bin)
      }
      
      # rarefy distacne data
      rare <- pblapply(nrar, function (n_samp) {
        reps <- lapply(seq_len(500), function (nrep) {
          # sample taxa
          tx <- sample(bin, size = n_samp, replace = FALSE)
  
          # create reduced matrices from dist_data[[n]]
          tdat <- upperTriangle(run$ddat[tx, tx])
          cdat <- upperTriangle(dist_data[[5]][tx, tx])
  
          # calculate pairwise distance metrics
          mpd <- mean(tdat, na.rm = TRUE)
          wpd <- sum(tdat * cdat, na.rm = TRUE) / sum(cdat, na.rm = TRUE)
          
          # return data frame
          data.frame(mpd = mpd,
                     wpd = wpd)
        }) %>%
          do.call(rbind, .)
        
        # sort values for disparity rarefactions
        smpd <- sort(reps$mpd)
        swpd <- sort(reps$wpd)
  
        # return data frames of rarefied values and confidence intervals
        list(mpd = data.frame(n = n_samp,
                   pd_mean = mean(smpd),
                   pd_low  = smpd[length(smpd) * 0.025],
                   pd_upp  = smpd[length(smpd) * 0.975]),
             wpd = data.frame(n = n_samp,
                   pd_mean = mean(swpd),
                   pd_low  = swpd[length(swpd) * 0.025],
                   pd_upp  = swpd[length(swpd) * 0.975]))
      })
      
      # return lists of values
      list(mpd = lapply(rare, "[[", "mpd") %>% do.call(rbind, .),
           wpd = lapply(rare, "[[", "wpd") %>% do.call(rbind, .))
    })
    
    # return lists of values
    list(list(rare = lapply(brare, "[[", "mpd"),
              dist = run$dist,
              bins = run$bins,
              disp = "mpd"),
         list(rare = lapply(brare, "[[", "wpd"),
              dist = run$dist,
              bins = run$bins,
              disp = "wpd"))
  })
}

scree <- function (mat) {
  # Calculates the variance contribution of each column of a matrix as a
  # percentage of the total, allowing production of a scree plot later on.
  #
  # Args:
  #   mat: the matrix on which to calculate, with variables (i.e. PC) as columns.
  #
  # Returns:
  #   A vector of percentage variance for each axis; if using PC(o) data then
  #   this will be in descending order.
  # get per-column variance
  pc_var <- apply(mat, 2, var)

  # convert to percentage of total
  pc_var <- pc_var / sum(pc_var) * 100

  # return
  return(pc_var)
}

subsetBootstrap <- function (mat, binning) {
  # Subsets and creates bootstrap replicates for a PC matrix using the functions
  # in dispRity. Requires packages `dispRity` and `magrittr`. `custom.subsets`
  # is used in preference to `chrono.subsets` as this allows prescription of
  # taxon occurrences across bin boundaries: taxa with FAD/LAD at a bin boundary
  # will occur in both bins, even if only intended for one.
  #
  # Args:
  #   mat:     the matrix to subset and bootstrap upon with variables (i.e. PC)
  #            as columns.
  #   binning: a list of bins each with the names of taxon occurrences matching
  #            the row names of `mat`; taxon names in `binning` but not in `mat`
  #            must be removed.
  #
  # Returns:
  #   A dispRity object with subsets, bootstrap and rarefaction replicates for
  #   use in subsequent disparity analyses.
  # subsetting and bootstrapping
  boot <- custom.subsets(mat,
                         group = binning) %>%
          boot.matrix(bootstraps = 500,
                      rarefaction = TRUE)

  # return data
  return(boot)
}

upperTriangle <- function (x, diag = FALSE, byrow = FALSE) {
  # Grabs only the upper triangle from a symmetrical matrix.
  # Taken from package `gdata`.
  #
  # Args:
  #   x: a symmetrical matrix
  #   diag: logical, include the diagonal values
  #   byrow: logical, return elements column- or row-wise
  #
  # Returns:
  #   A matrix of the upper traingle of matrix x.
  if (byrow)
    t(x)[rev(upper.tri(x, diag = diag))]
  else x[upper.tri(x, diag = diag)]
}
