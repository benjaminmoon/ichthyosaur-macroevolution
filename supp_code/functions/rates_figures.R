#!/usr/bin/env Rscript

# Rates figure plotting functions.
# From Moon & Stubbs "Early high rates and disparity in the evolution of
# ichthyosaurs. Communications Biology 

library(colorspace)

# Contents:
#   - fig_diversity
#   - fig_spaghettiHedman

fig_diversity <- function (diversity) {
  # Plot taxic and phylogenetic diversity from output put of divTaxPhy.
  #
  # Args:
  #   diversity: output pf divTaxPhy: a data.frame with columns "tex_div",
  #   "phy_meddiv", "phy_maxdiv", "phy_mindiv".
  #
  # Returns:
  #   Plots of taxic and phylogenetic diversity for both binning schemes.
  # plot for epochs
  plotDisparity(bins$epoch,
                diversity$epoch$tax_div,
                diversity$epoch[, c("phy_meddiv", "phy_maxdiv", "phy_mindiv")],
                xlim = range(bins))
  text(par()$usr[2], y = par()$usr[4], "a",
       adj = c(1.2, 1.2), font = 2, cex = 1.2)
  # plot for 10Ma bins
  plotDisparity(bins$tenMa_bJ,
                diversity$tenMa_bJ$tax_div,
                diversity$tenMa_bJ[, c("phy_meddiv", "phy_maxdiv", "phy_mindiv")],
                xlim = range(bins))
  text(par()$usr[2], y = par()$usr[4], "b",
       adj = c(1.2, 1.2), font = 2, cex = 1.2)

  # draw x-axis with periods
  xAxisPeriods()

  # plot axis labels
  mtext("Age (Ma)", side = 1, line = 3, outer = TRUE,
        font = 2, cex = 1.2)
  mtext("Diversity", side = 2, line = 2, outer = TRUE,
        font = 2, cex = 1.2)

  # add a legend – Tom Stubbs – in an overlay
  disparityLegend("bottomleft",
                  "Taxic diversity",
                  "Median phylogenetic diversity",
                  cex = 0.8)
}

fig_spaghettiHedman <- function (runs) {
  # Plot spaghetti plots of discrete character rates on trees from the Hedman
  # scaling method.
  #
  # Args:
  #   runs: the set of outputs from DiscreteCharacterRates to use/include in the
  #         spaghetti plot.
  #
  # Returns:
  #   A speghetti plot of the requested DiscreteCharacterRate runs with both
  #   epoch- and 10 Ma-bins.
  # plot epoch bins
  list.find(runs,
            bins == "epoch",
            n = length(discrete_rates)) %>%
  base_spaghettiPlot(xlim = range(bins))
  text(par()$usr[2], y = par()$usr[4], "a",
       adj = c(1.2, 1.2), font = 2, cex = 1.2)
  # plot for 10Ma bins
  list.find(runs,
            bins == "tenMa_bJ",
            n = length(discrete_rates)) %>%
  base_spaghettiPlot(xlim = range(bins))
  text(par()$usr[2], y = par()$usr[4], "b",
       adj = c(1.2, 1.2), font = 2, cex = 1.2)

  # draw x-axis with periods
  xAxisPeriods()

  # plot axis labels
  mtext("Age (Ma)", side = 1, line = 3, outer = TRUE,
        font = 2, cex = 1.2)
  mtext(expression("Rate of character changes" ~ (Ma ^ -1)),
        side = 2, line = 1.5, outer = TRUE,
        font = 2, cex = 1.2)

  # add legend
  base_spaghettiLegend("bottomleft")
}
