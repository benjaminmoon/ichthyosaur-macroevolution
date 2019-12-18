#!/usr/bin/env Rscript

# Useful plotting functions for time series.
# From Moon & Stubbs "Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

library(colorspace)

# Contents:
#   - diverge_pal
#   - base_spaghettiLegend
#   - base_spaghettiPlot
#   - boxBins
#   - disparityLegend
#   - plotDisparity
#   - plotConvexHull
#   - plotXYMorphospace
#   - xAxisPeriods
#   - xyLegend


diverge_pal <- diverging_hcl(7, h = c(80, 265))[c(2, 7)]

base_spaghettiLegend <- function(location = "bottomleft") {
  # Plots a legend for `base_spaghettiPlot` in the bottom left (by default) of
  # the figure area. Can accommodate multi-panel plots. Requires `magrittr`.
  library(magrittr)

  # save old parameters
  oldPar <- par()[c("fig", "oma", "mar")]

  # text for legend
  ltxt <- c("Rate of character changes",
            "Mean rate of character changes",
            "Significantly high",
            "Significantly low")

  # set plot to margin-less
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

  # plot legend with lines
  legend(location,
         ncol = 2,
         legend = ltxt,
         text.col = "#00000000",
         col = c("#00000030", "#999999", "#00000030", "#00000030"),
         lwd = c(1, 5, 1, 1),
         lty = c(1, 1, 1, 1),
         pch = c(NA, NA, NA, NA),
         bty = "n", xpd = TRUE, cex = 0.8)
  # plot legend with points
  legend(location,
         ncol = 2,
         legend = ltxt,
         col = c("#000000", NA, diverge_pal[1], diverge_pal[2]),
         lwd = c(1, 5, 1, 1),
         lty = c(0, 0, 0, 0),
         pch = c(18, NA, 17, 19),
         bty = "n", xpd = TRUE, cex = 0.8)

  # reinstate original `par` settings
  on.exit(par(oldPar))
}

base_spaghettiPlot <- function (data, xlim = NULL) {
  # Produce a spaghetti plot using base R plotting functions. Takes the output
  # of Claddis::DiscreteCharacterRate. Reimplemented from the code of Close
  # _et al._ (2015, _Curr. Biol.,_ doi:10.1016/j.cub.2015.06.047). Supplementary
  # information found at http://doi.org/10.5061/dryad.760sc.
  #
  # Args:
  #   data:  list of outputs from Claddis:DiscreteCharaterRate. The matrix
  #          per.bin.rates is used to plot from.
  #   xrang: the x-axis range, for plotting multiple on one x-axis.
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
                             # return
                             return(m)
                           }),
                    mean = sapply(rates, function (run) {
                             # separate in.rate
                             m <- run[["data"]][, "in.rate"]
                             # return
                             return(m)
                           }) %>%
                           cbind %>%
                           rowMeans %>%
                           cbind(midpoints  = rates[[1]][["data"]][, "midpoints"],
                                 mean_rates = .) %>%
                           as.data.frame)

  # set plot area
  plot(plot_vals$xval, plot_vals$yval, type = "n", axes = FALSE,
       xlim = plot_vals$xmax, ylim = range(0, plot_vals$ymax), ann = FALSE)
  boxBins(rates[[1]][["data"]][, c("from", "to")])
  # plot spaghetti lines and points
  lapply(rates, function (run) {
    lines(run$data$midpoints, run$data$in.rate, col = "#00000030")
  }) %>% invisible
  lapply(rates, function (run) {
    points(run$ns$midpoints, run$ns$in.rate, pch = 18, col = "#000000")
    points(run$hi$midpoints, run$hi$in.rate, pch = 17, col = diverge_pal[1])
    points(run$lo$midpoints, run$lo$in.rate, pch = 19, col = diverge_pal[2])
  }) %>% invisible
  # plot per-bin mean line
  lines(plot_vals$mean[, "midpoints"], plot_vals$mean[, "mean_rates"],
        col = "#999999", lwd = 5)

  # add y axis
  axis(side = 2)

  # add box around plot
  box(which = "plot")
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
  for (bin in seq(1, nrow(bins), 2)) {
    rect(bins[bin, 1], par()$usr[4], bins[bin, 2], par()$usr[3],
         density = NA, col = "grey80", xpd = FALSE)
  }
  ## Add an extra bin if the above is even to clarify where the last bin ends
  if (nrow(bins) %% 2 == 0) {
    rect(min(bins), par()$usr[4], par()$usr[2], par()$usr[3],
         density = NA, col = "#cccccc", xpd = FALSE)
  }
}

disparityLegend <- function (location = "bottomleft", namA, namB = NULL, cex = 1.0) {
  # Plots a legend for the function `plotDisparity` in the bottom left (by
  # default) of the plot as an overlay. Incorporates single or dual time series.
  #
  # Args:
  #   location:   the location of the legend, bottom left by default.
  #   namA, namB: the names of the series plotted.
  #   cex:        character expansion for the legend.
  #
  # Returns:
  #   A legend with up to two series named and coloured, matching that of
  # `plotDisparity` (colours from package `colorspace`).
  # save old parameters
  oldPar <- par()[c("fig", "oma", "mar", "mfrow", "mfcol", "cex")]

  # start new overlay, set plot to margin-less
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

  # plot the legend
  legend(location,
         legend = c(namA, namB),
         ncol = 1,
         pch = c(16, 18),
         col = c(diverge_pal[1], diverge_pal[2]),
         bty = "n", xpd = TRUE, cex = cex)

  # reinstate original `par` settings
  on.exit(par(oldPar))
}

plotDisparity <- function (binning, datA, datB = NULL, xlim = NULL, two_axes = FALSE) {
  # Plot two binned time series with error bars.
  #
  # Args:
  #   binning:    data.frame with two columns of bin start and end dates.
  #   datA, datB: data frames with one or three columns: disparity summary
  #               metric (first) with optionally lower and upper range/error
  #               intervals.
  #   xlim:       vector of values on which to calculate x-axis limits; use to
  #               plot several panels with same age range.
  #   two_axes:   logical, whether to overlay plots on separate y-axes
  #
  # Returns:
  #   Plots grey boxes for alternating bins, time series lines, error bars.
  #   get bin means.
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
       ann  = FALSE, axes = FALSE, type = "n")
  boxBins(binning)

  # plot first time series and error bars in datA
  if (ncol(datA) == 3) {
    arrows(plot_x + x_shift, datA[, 2], plot_x + x_shift, datA[, 3],
           length = 0.02, angle = 90, code = 3, col = diverge_pal[1])
    }
  points(plot_x + x_shift, datA[, 1], cex = 1.2, lwd = 1.2,
         type = "o", pch = 16, col = diverge_pal[1])
  axis(side = 2)

  # plot second time series and error bars in datB, if present
  if (!is.null(datB)) {
    # plot with separate y-axis
    if (two_axes == TRUE) {
      par(new = TRUE)
      plot(plot_x - x_shift, datB[, 1],
           xlim = rev(x_range),
           ylim = yB_range,
           ann  = FALSE, axes = FALSE, type = "n")
      axis(side = 4)
    }
    if (ncol(datB) == 3) {
      arrows(plot_x - x_shift, datB[, 2], plot_x - x_shift, datB[, 3],
             length = 0.02, angle = 90, code = 3, col = diverge_pal[2])
    }
    points(plot_x - x_shift, datB[, 1], cex = 1.2, lwd = 1.2,
           type = "o", pch = 18, col = diverge_pal[2])

  }

  # surround plot with a box
  box(which = "plot")
}

plotConvexHull <- function (xcoord, ycoord, lcol){
  # Plots a 2D minimal convex hull around the given set of points from x- and
  # y-coordinates, in the desired colour.
  #
  # Args:
  #   xcoord, ycoord: matching x- and y-coordinates of the points to surround.
  #   lcol:           the colour of the line to draw.
  #
  # Returns:
  #   On a plot returns the minimal convex hull for the specified points with
  #   the specified line colour.
  # computes convex hull from points giving the points demarcating the convex
  # hull
  hpts <- chull(x = xcoord, y = ycoord)

  # adds first point to join the hull
  hpts <- c(hpts, hpts[1])

  # plots the coordinates of the hull-demarcating points and joins by lines
  lines(xcoord[hpts], ycoord[hpts], col = lcol)
}

plotXYMorphospace <- function (pcdat, bindat, bincol, binpch = 19, axes = TRUE) {
  # Plots an x-y scatter plot, divided into bins and with minimal convex hulls.
  # Can be integrated into multiple plots. Points are printed after convex hulls
  # so they overlie all lines. No axes are labelled.
  #
  # Args:
  #   pcdat:  a matrix or data.frame with row names and coordinates in three
  #           columns; only the first three columns are used if there are more.
  #   bindat: a list of elements, each element containing row names for each bin.
  #   bincol: a vector/list of colours matching the length of bindat.
  #   binpch: plotting character for each bin, default is all 19 (large circle).
  #   axes:   logical, print axes or not?, useful for inclusion in multi-panel
  #           figures.
  #
  # Returns:
  #   A figure with two plots arranged vertically, sharing the same x-axis
  #   (PCo 1), and points coloured and plotted according to their bin.
  # begin plot
  if (axes == FALSE) {
    # without axes
    plot(pcdat[, 1],
         pcdat[, 2],
         type = "n",
         ann = FALSE, xaxt = "n", yaxt = "n")
  } else {
    # with axes
    plot(pcdat[, 1],
         pcdat[, 2],
         type = "n")
  }

  # plot convex hulls for each bin
  for (bin in seq_along(bindat)) {
    # subset points to plot
    dat <- pcdat[row.names(pcdat) %in% bindat[[bin]], ]

    # plot convex hull
    plotConvexHull(xcoord = dat[, 1],
                   ycoord = dat[, 2],
                   lcol   = bincol[[bin]])
  }

  # plot points for each bin
  for (bin in seq_along(bindat)) {
    # subset points to plot
    dat <- pcdat[row.names(pcdat) %in% bindat[[bin]], ]

    # plot points
    points(dat[, 1],
           dat[, 2],
           col = bincol[[bin]], pch = binpch[bin])
  }
}

xAxisPeriods <- function (cex = 1.0) {
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
  axis(1, pos = axis.bottom)
  # plot period boxes
  # Triassic
  rect(min(par()$usr[1], 252.17), axis.bottom, 201.3, axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(min(par()$usr[1], 252.17), 201.3)),
       mean(c(axis.bottom, axis.top)),
       "TRIASSIC",
       adj = c(0.5, 0.5), cex = cex, xpd = NA)
  # Jurassic
  rect(201.3, axis.bottom, 145, axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(201.3, 145)),
       mean(c(axis.bottom, axis.top)),
       "JURASSIC",
       adj = c(0.5, 0.5), cex = cex, xpd = NA)
  # Cretaceous
  rect(145, axis.bottom, max(par()$usr[2], 66.0), axis.top,
       col = "white", border = TRUE, xpd = NA)
  text(mean(c(145, max(par()$usr[2], 66.0))),
       mean(c(axis.bottom, axis.top)),
       "CRETACEOUS",
       adj = c(0.5, 0.5), cex = cex, xpd = NA)
  # tidy up box around
  rect(par()$usr[1], axis.bottom, par()$usr[2], axis.top,
       col = NA, border = TRUE, xpd = NA)
}

xyLegend <- function (location = "bottom", binnam, bincol, binpch = 19, colnum = length(binnam), cex = 1.0) {
  # Plots a legend for the function `plotXYMorphospace` in the bottom (by
  # default) of the plot as an overlay.
  #
  # Args:
  #   location: the location of the legend, bottom left by default.
  #   binnam:   vector, the names of the bins plotted.
  #   bincol:   vector, the colours of the plotted bins
  #   binpch:   vector, the printing characters of the plotted bins, defaults to
  #             19 (large circle).
  #   colnum:   number of columns.
  #
  # Returns:
  #   A legend with up to two series named and coloured, matching that of
  # `plotDisparity` (colours from package `colorspace`).
  # save old parameters
  oldPar <- par()[c("fig", "oma", "mar")]

  # start new overlay, set plot to margin-less
  par(fig = c(0, 1, 0, 1),
      oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0),
      new = TRUE)
  plot(0, 0,
       type = "n",
       bty  = "n",
       xaxt = "n",
       yaxt = "n")

  # plot the legend
  legend(location,
         legend = binnam,
         ncol   = colnum,
         pch    = binpch,
         col    = bincol,
         bty = "n", xpd = TRUE, cex = cex)

  # reinstate original `par` settings
  on.exit(par(oldPar))
}
