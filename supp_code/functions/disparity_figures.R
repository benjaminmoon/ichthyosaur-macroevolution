#!/usr/bin/env Rscript

# Disparity figure plotting functions.
# From Moon & Stubbs <Title>, <doi>

library(colorspace)

source("functions/plotting_functions.R")

# Contents:
#   - fig_disparityAll
#   - fig_maxPD
#   - fig_pdAll
#   - fig_rarefactionCurves
#   - fig_scree
#   - fig_xyMorphospace
#   - fig_xyMorphospaceAll

fig_disparityAll <- function () {
  # Plot disparity-time curves for all analyses
  # top level: metrics (intended to plot one-per-page)
  for (disp in names(plot_names$disp_metrics)) {
    for (distance in plot_names$dist) {
      # circle through binning schemes (in single columns)
      for (scheme in names(bins)) {
        # get data
        dat <- list.find(disparity,
                         metric == disp & dist == distance & bins == scheme,
                         n = length(disparity))

        # plot disparity curves
        plotDisparity(bins[[scheme]],
                      list.find(dat,
                                corr == "uncorrected")[[1]]$summ_test[,
                                                                       c("bs.mean",
                                                                         "2.5%",
                                                                         "97.5%")],
                      list.find(dat,
                                corr == "Caillez-corrected")[[1]]$summ_test[,
                                                                             c("bs.mean",
                                                                               "2.5%",
                                                                               "97.5%")],
                      two_axes = TRUE, xlim = range(bins))

        # plot distance (column) title
        if (par()$fig[4] == 1) mtext(distance, side = 3, font = 2, cex = 0.8)
      }
      # print x-axis
      xAxisPeriods(cex = 0.8)
    }

    # print titles
    mtext(paste0("Ichthyosaur disparity (mean ",
                 plot_names$disp_metrics[[disp]]$name,
                 ") through time"),
          outer = TRUE, line = 2, font = 2)
    mtext("Age (Ma)", side = 1, font = 2, outer = TRUE, line = 4, cex = 0.9)
    mtext(paste("Mean",
                plot_names$disp_metrics[[disp]]$name,
                "(uncorrected PCo)"),
          side = 2, font = 2, outer = TRUE, cex = 0.9)
    mtext(paste("Mean",
                plot_names$disp_metrics[[disp]]$name,
                "(corrected PCo)"),
          side = 4, font = 2, outer = TRUE, cex = 0.9)

    # plot legend
    disparityLegend(location = "bottomleft",
                    paste(plot_names$disp_metrics[[disp]]$name,
                          "(uncorrected PCo)"),
                    paste(plot_names$disp_metrics[[disp]]$name,
                          "(Caillez-corrected PCo)"),
                    cex = 0.8)
    # prepare for next page
    par(mfcol = c(2, 4))
  }
}

fig_maxPD <- function (metric) {
  # Plot the disparity as pairwise distances of the MAX distance matrix.
  #
  # Args:
  #   metric: the metric to be used: a list combining data.frames of MPD and WMPD
  #           with three columns (mean, lower and upper) each, for two binning schemes.
  #
  # Returns:
  #   Two disparity-time plots in a column, one for each binning scheme, and
  #   each plotting MPD and WMPD with error bars.
  # plot epoch time series
  plotDisparity(bins$epochs, metric$epochs$MPD, metric$epochs$WMPD,
                xlim = range(bins))

  # label for metric
  text(par()$usr[2],
       par()$usr[4],
       "a",
       font = 2, adj = c(1.2, 1.2))

  # plot equal time series
  plotDisparity(bins$tenMa_bJ, metric$tenMa_bJ$MPD, metric$tenMa_bJ$WMPD,
                xlim = range(bins))

  # label for metric
  text(par()$usr[2],
       par()$usr[4],
       "b",
       font = 2, adj = c(1.2, 1.2))

  # plot x-axis
  xAxisPeriods()

  # plot axis labels
  mtext("Age (Ma)", side = 1, line = 3, outer = TRUE, cex = 1.2, font = 2)
  mtext("Pairwise maximum observed rescaled distance",
        side = 2, line = 2, outer = TRUE, cex = 1.2, font = 2)

  # add a legend: Tom Stubbs
  disparityLegend("bottomleft",
                  "Pairwise maximum observed rescaled distance",
                  "Weighted pairwise maximum observed rescaled distance",
                  cex = 0.8)
}

fig_pdAll <- function () {
  # Plot all pairwise distances disparity plots
  for (distNam in 1:4) {
    # get required data
    metric <- pd_disparity[[distNam]]

    # plot column
    plotDisparity(bins$epochs, metric$epochs$MPD, metric$epochs$WMPD,
                  xlim = range(bins))

    # label
    mtext(paste0(plot_names$dist[distNam]),
          font = 2, side = 3, cex = 0.8)

    # plot equal bins
    plotDisparity(bins$tenMa_bJ, metric$tenMa_bJ$MPD, metric$tenMa_bJ$WMPD,
                  xlim = range(bins))

    # plot x-axis labels
    xAxisPeriods()
  }

  # plot for character comparisons; weighting doesn't make sense here
  metric <- pd_disparity$comp.char.matrix

  # plot epochs
  plotDisparity(bins$epochs, metric$epochs$MPD,
                xlim = range(bins))

  # label
    mtext("CHAR",
          font = 2, side = 3, cex = 0.8)

  # plot equal bins
  plotDisparity(bins$tenMa_bJ, metric$tenMa_bJ$MPD,
                xlim = range(bins))

  # plot x-axis labels
  xAxisPeriods()

  # print titles
  mtext(paste0("Ichthyosaur disparity (mean pairwise dissimilarity) through time"),
        outer = TRUE, line = 2, font = 2, side = 3)
  mtext("Age (Ma)", side = 1, line = 3.54, outer = TRUE,
        adj = c(0.5, 0.5), cex = 1)
  mtext("Pairwise distance", side = 2, line = 0.5, outer = TRUE,
        adj = c(0.5, 0.5), cex = 1)

  # legend: Stubbsy
  disparityLegend("bottomleft",
                  "Pairwise distance",
                  "Weighted pairwise distance")
}

fig_rarefactionCurves <- function (disparity) {
  # Plot rarefaction curves for each bin; new page for each run.
  for (run in disparity) {
    # get disparity data
    sum <- run$disp_summ

    # get subsets
    lev <- unique(run$disp_summ$subsets)

    # prepare plot
    par(mfrow = c(ceiling(length(lev) / 3), 3))

    # plot per level (bin)
    for (level in lev) {
      # get rows to plot
      set <- which(sum$subsets == level)
      # get y limits
      yl <- range(sum[set, c("bs.mean", "2.5%", "97.5%")])
      # where no taxa present
      if (is.na(yl[1])) {
        yl <- range(0, 1)
      }
      # plot count against disparity metric
      plot(sum$n[set],
           sum$bs.mean[set],
           type = "n", ann = FALSE,
           xlim = range(0, max(sum$n)),
           ylim = yl)

      # add confidence interval polygon
      polygon(c(sum$n[set], rev(sum$n[set])),
              c(sum$`2.5%`[set], rev(sum$`97.5%`[set])),
              border = NA,
              col = "grey70")

      # plot mean disparity points
      points(sum$n[set],
             sum$bs.mean[set],
             pch = 19, type = "o")

      # text "n="
      text(par()$usr[2], par()$usr[4],
           paste("n =", max(sum$n[set])),
           adj = c(1.2, 1.2))

      # plot title; if epoch bins give epoch name
      if (run$bins == "epochs") {
        title(main = plot_names$epoch_names[[which(level == lev)]],
              line = 0.5)
      } else {
        title(main = level, line = 0.5)
      }
    }

    # add titles
    mtext(paste("Rarefaction curves: mean",
                plot_names$disp_metrics[[run$metric]]$name,
                "of",
                run$corr, run$dist,
                "distance matrix in",
                plot_names$bins[run$bins],
                "bins"),
          side = 3, outer = TRUE, font = 2)
    mtext("Taxon count", side = 1, outer = TRUE, font = 2, line = 1)
    mtext(paste("Mean", plot_names$disp_metrics[[run$metric]]$name),
          side = 2, outer = TRUE, font = 2, line = 2)
  }
}

fig_scree <- function () {
  # Scree plots of PCo analyses on distance matrices and negative eigenvalue
  # corrections.
  # max data length
  max_axes <- lapply(scree_data, "[[", 1) %>%
              lengths %>%
              max

  # colours for lines
  pal <- rainbow_hcl(4)

  # plot scree data
  # repeat for each item in `scree_data`
  for (mat in seq_along(scree_data)) {
    # get y-axis limits per plotting row
    if (mat < 5) {
      yrang <- range(0, max(unlist(scree_data[1:4])))
    } else {
      yrang <- range(0, max(unlist(scree_data[5:8])))
    }

    # plotting values for columns or rows
    col <- mat %% 5
    row <- (mat %/% 5) + 1

    # plot scree values
    plot(scree_data[[mat]]$scree, type = "l", lty = 1, lwd = 2,
         xlim = range(0, max_axes), ylim = yrang,
         col = rep(pal, 2)[mat], xaxt = "n", yaxt = "n")
    # plot y-axis to left
    if (par()$fig[1] < 0.25) axis(2)

    # plot PCo correlations
    par(new = TRUE)
    plot(scree_data[[mat]]$pc_cor, pch = 16,
         xlim = range(0, max_axes), ylim = range(0, 1),
         col = rep(pal, 2)[mat], xaxt = "n", yaxt = "n")

    # plot axes
    if (par()$fig[2] > 0.75) axis(4)
    if (par()$fig[3] < 0.5) axis(1)

    # titles to plots
    # distance matrix
    if (par()$fig[4] > 0.5) {
      mtext(plot_names$dist[col],
            side = 3, line = 0.5, adj = c(0.5, 0.5), cex = 0.9)
    }
    # Caillez correction
    if (par()$fig[1] < 0.25) {
      mtext(plot_names$corr[row],
            side = 2, line = 1.5, adj = c(0.5, 0.5), cex = 0.9)
    }
  }

  # figure and axis titles
  mtext("Scree plots of PCo variance", outer = TRUE, line = 2, side = 3, font = 2)
  mtext("PC axis", outer = TRUE, line = 2, side = 1, font = 2)
  mtext("Percentage variance", outer = TRUE, line = 3, side = 2, font = 2)
  mtext("Correlation between distance matrix and n PCo axes",
        outer = TRUE, line = 2, side = 4, font = 2)

  # add legend
  # set plot to margin-less
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottomleft",
         ncol = 2,
         legend = c("Variance", "Correlation"),
         lty = c(1, NA),
         lwd = 2,
         pch = c(NA, 16),
         bty = "n", xpd = TRUE)
}

fig_xyMorphospace <- function () {
  # Plot the x-y morphospace for PCo 1–3 of the uncorrected PCo analysis on MAX
  # distance matrix.
  # plot PCo 1–3 scatter plots
  for (pc_axis in PC$y) {
    # plot xy morphospace
    plotXYMorphospace(pco_data[[pc_matrix]]$pco[, c(1, pc_axis)],
                      bin_data$epoch$bin_data[-9],
                      epoch_col,
                      epoch_pch, axes = FALSE)

    # plot y-axis
    axis(2)
    mtext(paste0("PCo ", pc_axis, " (",
                 signif(pco_data[[pc_matrix]]$scree[pc_axis], 3),
                 "% variance)"),
          side = 2, line = 1.5, font = 2)

    # plot metric above (only if top plot)
    text(par()$usr[2], par()$usr[4],
         c("MAX", "a", "b")[pc_axis],
         font = 2, adj = c(1.2, 1.2))
  }

  # plot x-axis
  axis(1)
  mtext(paste0("PCo ", PC$x, " (",
               signif(pco_data[[pc_matrix]]$scree[PC$x], 3),
               "% variance)"),
        side = 1, line = 2, font = 2)

  # add legend: Tom Stubbs
  xyLegend(binnam = plot_names$epoch_names,
           bincol = unlist(epoch_col),
           binpch = epoch_pch,
           colnum = 3,
           cex = 0.7)
}

fig_xyMorphospaceAll <- function () {
  # Plot x-y morphospace for PCo 1–3 of all PCo analyses.
  # plot x-y scatter plots for PCo 1–3
  for (pc_matrix in seq_along(pco_data)) {
    # plot xy and xz for PCo axes 1–3
    for (pc_axis in PC$y) {
      # plot xy morphospace
      plotXYMorphospace(pco_data[[pc_matrix]]$pco[, c(1, pc_axis)],
                        bin_data$epoch$bin_data[-9],
                        epoch_col,
                        epoch_pch, axes = FALSE)

      # plot y-axis
      axis(2)
      mtext(paste0("PCo ", pc_axis, " (",
                   signif(pco_data[[pc_matrix]]$scree[pc_axis], 3),
                   "% variance)"),
            side = 2, line = 1.5, cex = 0.7, font = 2)

      # plot metric above (only if top plot)
      if (par()$fig[4] > 0.75) {
        mtext(rep(plot_names$dist, 2)[pc_matrix], side = 3, font = 2, cex = 0.9)
      }
    }

    # plot x-axis
    axis(1)
    mtext(paste0("PCo ", PC$x, " (",
                 signif(pco_data[[pc_matrix]]$scree[PC$x], 3),
                 "% variance)"),
          side = 1, line = 2, cex = 0.7, font = 2)

    # add title
    if (pc_matrix == 4) {
      mtext("PCo 1-3 scatter plots without negative eigenvalue correction",
            outer = TRUE, line = 2, font = 2)
    }
    if (pc_matrix == 8) {
      mtext("PCo 1-3 scatter plots with negative eigenvalue correction",
            outer = TRUE, line = 2, font = 2)
    }

    # add legend when figure is full
    if (par()$fig[2] == 1 & par()$fig[4] == 0.5) {
    xyLegend(binnam = plot_names$epoch_names,
             bincol = unlist(epoch_col),
             binpch = epoch_pch)
    par(mfcol = c(2, 4))
    }
  }
}
