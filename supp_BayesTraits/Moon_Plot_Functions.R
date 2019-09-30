library(colorspace)
library(magrittr)
library(rlist)

diverge_pal <- diverge_hcl(7, h = c(340, 128))[c(2, 7, 3)]

base_spaghettiPlot <- function (data, col = diverge_pal, xlim = NULL, cex_ax = 1.0) {
  # Produce a spaghetti plot using base R plotting functions. Takes the output
  # of Claddis::DiscreteCharacterRate. Reimplemented from the code of Close
  # _et al._ (2015, _Curr. Biol.,_ doi:10.1016/j.cub.2015.06.047). Supplementary
  # information found at http://doi.org/10.5061/dryad.760sc.
  #
  # Args:
  #   data:  list of outputs from Claddis:DiscreteCharaterRate. The matrix
  #          per.bin.rates is used to plot from.
  #   col:   vector of three colours, by default taken from diverge_pal (defined
  #          above).
  #   xrang: the x-axis range, for plotting multiple on one x-axis.
  #   cex:   character expansion.
  #
  # Returns:
  #   A plot of per-bin rates for each tree (spaghetti plot) with coloudatA
  #   points for significatly high/low per-bin rates (datA/datB), and a thick
  #   grey mean line. Can be incorporated into multi-panel plots. A legend is
  #   not included, but can be added separately.
  # get per.bin.rates matrix
  rates <- lapply(data, function (run) {
    # get data and convert to data frame
    df <- run[["per.bin.rates"]] %>%
          as.data.frame

    # calculate bin midpoints
    df$midpoints <- rowMeans(df[, c("from", "to")])

    # get values for significantly high/low rates and non-significant
    hi <- df[which(df$ml.signif.hi == 1), c("midpoints", "in.rate")]
    lo <- df[which(df$ml.signif.lo == 1), c("midpoints", "in.rate")]
    ns <- df[which(df$ml.signif == 0), c("midpoints", "in.rate")]

    # gather into a list
    return(list(data = df,
                hi   = hi,
                lo   = lo,
                ns   = ns))
  })

  # x-axis range
  if (is.null(xlim)) {
    xlim <- range(rates[[1]][["data"]][, c("from", "to")])
  }

  # maximum rates and mean rates
  plot_vals <- list(xval = rates[[1]][["data"]][, "midpoints"],
                    yval = rates[[1]][["data"]][, "in.rate"],
                    xmax = rev(xlim),
                    ymax = sapply(rates, function (run) {
                                   # maximum per-bin rates
                                   m <- run[["data"]][, "in.rate"] %>%
                                        max
                                  }),
                    mean = sapply(rates, function (run) {
                             m <- run[["data"]][, "in.rate"]
                                  }) %>%
                           cbind %>%
                           rowMeans %>%
                           cbind(midpoints  = rates[[1]][["data"]][, "midpoints"],
                                 mean_rates = .) %>%
                           as.data.frame)

  # set plot area
  plot(plot_vals$xval, plot_vals$yval, type = "n", axes = FALSE,
       xlim = plot_vals$xmax, ylim = range(0, plot_vals$ymax), ann = FALSE)
  # boxBins(rates[[1]][["data"]][, c("from", "to")])
  # plot spaghetti lines and points
  lapply(rates, function (run) {
    lines(run$data$midpoints, run$data$in.rate, col = "#00000030")
  }) %>% invisible
  lapply(rates, function (run) {
    points(run$ns$midpoints, run$ns$in.rate, pch = 18, col = "#000000")
    points(run$hi$midpoints, run$hi$in.rate, pch = 17, col = col[1])
    points(run$lo$midpoints, run$lo$in.rate, pch = 19, col = col[2])
  }) %>% invisible
  # plot per-bin mean line
  lines(plot_vals$mean[, "midpoints"], plot_vals$mean[, "mean_rates"],
        col = col[3], lwd = 5, lend = 0)

  # add y axis
  axis(side = 2, cex.axis = cex_ax)

  # add box around plot
  # box(which = "plot")
}

boxBins <- function (bins) {
  # Plots alternating grey boxes for a given set of bins.
  #
  # Args:
  #   bins: data.frame with two columns of start and end dates.
  #
  # Returns:
  #   Grey boxes within plots.
  # Get bin length for the data
  # for (bin in seq(1, nrow(bins), 2)) {
  #   rect(bins[bin, 1], par()$usr[4], bins[bin, 2], par()$usr[3],
  #        density = NA, col = "#cccccc", xpd = FALSE)
  # }
  rect(bins[, 1],
       par()$usr[3],
       bins[, 2],
       par()$usr[4],
       col = c("#cccccc", "#ffffff"),
       border = NA,
       xpd = FALSE)
  ## Add an extra bin if the above is even to clarify where the last bin ends
  if (nrow(bins) %% 2 == 0) {
    rect(bins[nrow(bins), 2],
         par()$usr[4],
         par()$usr[2],
         par()$usr[3],
         density = NA,
         col = "#cccccc",
         xpd = FALSE)
  }
}

plotDisparity <- function (binning, datA, datB = NULL, col = diverge_pal, xlim = NULL, two_axes = FALSE, cex_ax = 1.0) {
  # Plot two binned time series with error bars.
  #
  # Args:
  #   binning:    data.frame with two columns of bin start and end dates.
  #   datA, datB: data frames with one or three columns: disparity summary
  #               metric (first) with optionally lower and upper range/error
  #               intervals.
  #   col:        vector of two colours to use, defaults to diverge_pal (defined
  #               above).
  #   xlim:       vector of values on which to calculate x-axis limits; use to
  #               plot several panels with same age range.
  #   two_axes:   logical, whether to overlay plots on separate y-axes
  #
  # Returns:
  #   Plots grey boxes for alternating bins, time series lines, error bars.
  # get bin means.
  plot_x <- rowMeans(binning)
  datA <- as.data.frame(datA)

  # get x-axis limits
  if (is.null(xlim)) {
    x_range <- range(binning)
  } else {
    x_range <- range(xlim)
  }

  # get y-axis limits
  if (is.null(datB)) {
    y_range <- range(datA, na.rm = TRUE)
    x_shift <- 0
  } else if (two_axes == FALSE) {
    y_range <- range(datA, datB, na.rm = TRUE)
    x_shift <- 1.5
    datB <- as.data.frame(datB)
  } else {
    y_range <- range(datA, na.rm = TRUE)
    yB_range <- range(datB, na.rm = TRUE)
    x_shift <- 1.5
    datB <- as.data.frame(datB)
  }

  # empty plot + bins
  plot(plot_x, datA[, 1],
       xlim = rev(x_range),
       ylim = y_range,
       ann  = FALSE, axes = FALSE, type = "n", bty = "n")
  boxBins(binning)

  # plot first time series and error bars in datA
  if (ncol(datA) == 3) {
    arrows(plot_x + x_shift, datA[, 2], plot_x + x_shift, datA[, 3],
           length = 0.02, angle = 90, code = 3, col = col[1])
    }
  points(plot_x + x_shift, datA[, 1],
         type = "o", pch = 16, col = col[1])
  axis(side = 2, cex.axis = cex_ax)

  # plot second time series and error bars in datB, if present
  if (!is.null(datB)) {
    # plot with separate y-axis
    if (two_axes == TRUE) {
      par(new = TRUE)
      plot(plot_x - x_shift, datB[, 1],
           xlim = rev(x_range),
           ylim = yB_range,
           ann  = FALSE, axes = FALSE, type = "n", bty = "n")
      axis(side = 4, cex.axis = cex_ax)
    }
    if (ncol(datB) == 3) {
      arrows(plot_x - x_shift, datB[, 2], plot_x - x_shift, datB[, 3],
             length = 0.02, angle = 90, code = 3, col = col[2])
    }
    points(plot_x - x_shift, datB[, 1],
           type = "o", pch = 18, col = col[2])

  }
}

timeScalePhylo <- function (tree, ages, binning, xlim, cex_tip = 1.0, width = 1.0, label_offset, show_tip_label = TRUE) {
  # A modified version of strap::geoscalePhylo that can be plotted alongside
  # other plots and without the time scale.
  #
  # Args:
  #   tree:    a phylo object to plot.
  #   ages:    a two-column matrix or data.frame of taxon FAD and LAD, with row
  #            names matching the tip names in tree.
  #   binning: a two-column matrix or data frame of bin start an end dates.
  #   x_lim:    a vector of dates to show on the x-axis.
  #   cex_tip: character expansion of the tip labels.
  #   width:   edge width.
  #   label_offset: offset for tip labels.
  #
  # Returns:
  #   A phylogenetic tree plot with terminal ranges for each taxon and grey
  #   boxes to show bin extents.
  # get ages matches to tip.label
  if (missing(ages) == FALSE) {
    ranges <- TRUE
  } else {
    ranges <- FALSE
  }
  if (!missing(ages)) ages <- ages[tree$tip.label, ]

  # get xlim if not specified
  root_age <- tree$root.time
  if (!missing(xlim)) {
    xlim <- sort(root_age - xlim)
  } else if (ranges == TRUE && !missing(ages) && missing(xlim)) {
    xlim <- (root_age - min(ages)) + diff(range(ages)) * 0.05
  } else {
    xlim <- NULL
  }

  # rescale time scale for drawing bin-boxes
  binning <- binning[order(binning[, 1], decreasing = T), ]
  bins_rescaled <- root_age - binning

  # get taxon ranges and label offset
  first_la <- tree$root.time - dist.nodes(tree)[1, Ntip(tree) + 1]
  if (ranges == TRUE && missing(ages) == FALSE) {
    offset <- array(dim = length(tree$tip.label), data = 1)
    offset_correction <- diff(range(ages)) * 0.01
    taxon_ranges <- root_age - ages[, c("FAD", "LAD")]
    if (first_la != ages[1, "LAD"]) {
      if (!missing(label_offset)) {
        offset <- array(dim = length(ages[, "FAD"]),
                        data = (ages[, "FAD"] - ages[, "LAD"]) + label_offset)
      } else {
        offset <- array(dim = length(ages[, "FAD"]),
                        data = (ages[, "FAD"] - ages[, "LAD"]) + offset_correction)
      }
    }
  } else if (!missing(label_offset)) {
    offset <- label_offset
  } else {
    offset <- 1
  }

  # get phylogeny plot limits
  plot.phylo(tree,
             plot = FALSE,
             no.margin = TRUE,
             x.lim = xlim,
             direction = "rightwards",
             cex = cex_tip,
             show.tip.label = show_tip_label)

  # get figure position
  par_mfg <- par()$mfg
  par_fig <- par()$fig
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  # plot grey bars for bins
  boxBins(bins_rescaled)

  # new plot over old
  par(mfg = par_mfg, new = TRUE)

  # plot shown tree
  plot.phylo(tree,
             plot = TRUE,
             no.margin = TRUE,
             x.lim = xlim,
             direction = "rightwards",
             cex = cex_tip,
             show.tip.label = show_tip_label)

  # plot taxon ranges
  par(lend = 1)
  segments(taxon_ranges[, "FAD"],
           lastPP$yy[c(1:length(tree$tip.label))],
           taxon_ranges[, "LAD"],
           lastPP$yy[c(1:length(tree$tip.label))],
           col = "black",
           lwd = width * 2)

  # restore original par
  par(mfrow = par_mfg[3:4],
      mfg   = par_mfg,
      new = FALSE)
}

xAxisPeriods <- function (cex_title = 1.0, cex_ax = 1.0) {
  # Plots boxes with geological period labels along the x-axis.
  # Only for periods of the Mesozoic so far.
  #
  # Args:
  #   cex: character expansion of labels
  #
  # Returns:
  #   In plots, an x-axis with boxes for each geological period, labels, and
  #   time on the x-axis.
  # get x-axis extent
  axis.bottom <- par()$usr[3] - 0.1 * diff(par()$usr[3:4])
  axis.top <- par()$usr[3]
  # plot x-axis
  axis(1, pos = axis.bottom, cex.axis = cex_ax)
  # plot period boxes
  # Triassic
  rect(min(par()$usr[1], 252.17), axis.bottom, 201.3, axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(min(par()$usr[1], 252.17), 201.3)),
       mean(c(axis.bottom, axis.top)),
       "TRIASSIC",
       adj = c(0.5, 0.5), cex = cex_title, xpd = NA)
  # Jurassic
  rect(201.3, axis.bottom, 145, axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(201.3, 145)),
       mean(c(axis.bottom, axis.top)),
       "JURASSIC",
       adj = c(0.5, 0.5), cex = cex_title, xpd = NA)
  # Cretaceous
  rect(145, axis.bottom, max(par()$usr[2], 66.0), axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(145, max(par()$usr[2], 66.0))),
       mean(c(axis.bottom, axis.top)),
       "CRETACEOUS",
       adj = c(0.5, 0.5), cex = cex_title, xpd = NA)
  # tidy up box around
  rect(par()$usr[1], axis.bottom, par()$usr[2], axis.top,
       col = NA, border = TRUE, xpd = NA)
}