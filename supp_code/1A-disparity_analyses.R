#!/usr/bin/env Rscript

# 1A. Disparity analyses #
# ====================== #

# Script from Moon & Stubbs <Title>, <doi>

# required libraries
library(dispRity)
library(doParallel)
library(magrittr)
library(parallel)
library(pbapply)
library(rlist)

# required scripts
source("functions/disparity_functions.R")
source("functions/disparity_figures.R")
source("functions/strat_bin_colours.R")

# for parallel processing
# register parallel cores for disparity loop below
# change number of cores as available/desired
no_cores <- detectCores() - 1
clus <- makeCluster(no_cores, outfile = "")
registerDoParallel(clus)


# 1.1. Gather data #
# ---------------- #

# get data
dist_data <- Claddis::ReadMorphNexus("data/matrix.nex") %>%
             Claddis::MorphDistMatrixFast()

# get taxon ages
taxon_ages <- read.table("./data/ichthyosaur_occurrences.tsv",
                         header    = TRUE,
                         row.names = 1,
                         sep       = "\t")

# get bins
bins <- list(# bins corresponding to epochs
             epochs = read.table("./data/epoch_bins.tsv",
                                 header    = TRUE,
                                 row.names = 1,
                                 sep       = "\t"),
             # 10 Ma bins aligned to the base of the Jurassic (201.3 Ma)
             tenMa_bJ = data.frame(max_age = seq(251.3, 100, -10),
                                   min_age = seq(241.3, 90, -10),
                                   row.names = paste0("bin", seq(1, 16))))

# bin taxa and remove outgroup taxon _Hupehsuchus_nanchangensis_
outgroup <- "Hupehsuchus_nanchangensis"
bin_data <- lapply(bins, function(scheme) {
  data <- list()
  for (bin in seq_along(rownames(scheme))) {
    data[[bin]] <- rownames(taxon_ages)[which(taxon_ages$FAD > scheme[bin, "min_age"] &
                                              taxon_ages$LAD < scheme[bin, "max_age"] &
                                              !(rownames(taxon_ages) %in% outgroup))]
  }
  names(data) <- rownames(scheme)
  data$outgroup <- outgroup
  return(list(bins = scheme, bin_data = data))
})
saveRDS(bin_data, "output/bin_data.rds")

# names for pretty plotting
plot_names <- list(dist = c("RAW", "GED", "GOW", "MAX"),
                   corr = c("uncorrected", "Caillez-corrected"),
                   bins = c(epochs   = "epoch-length",
                            tenMa_bJ = "10 Ma"),
                   disp_metrics = list(sum_var  = list(func = c(sum, variances),
                                                       name = "sum of variances"),
                                       sum_rang = list(func = c(sum, ranges),
                                                       name = "sum of ranges"),
                                       centroid = list(func = c(centroids),
                                                       name = "centroid distance")),
                   epoch_names = c("Early Triassic",
                                   "Middle Triassic",
                                   "Late Triassic",
                                   "Early Jurassic",
                                   "Middle Jurassic",
                                   "Late Jurassic",
                                   "Early Cretaceous",
                                   "Late Cretaceous"))
# list of disparity metrics
clusterExport(clus, c("plot_names"))


# 1.2. Pairwise distances #
# ----------------------- #

# perform pairwise distance analyses
# creates a list: distance matrices, each containing a list: binning schemes,
# each containing a list: MPD and WMPD data frames
clusterExport(clus, c("bootstrapMPD", "bootstrapWMPD", "doPD", "upperTriangle"))
pd_disparity <- pblapply(dist_data, function(dissim) {
                  lapply(bin_data, function(scheme) {
                    doPD(scheme$bin_data[scheme$bin_data != scheme$bin_data$outgroup],
                         dissim,
                         dist_data$comp.char.matrix)
                  })
                })
saveRDS(pd_disparity, file = "output/pd_disparity.rds")

{
  # plot single distance metric
  metric <- pd_disparity$max.dist.matrix
  cairo_pdf("fig/fig1-pairwise_max.pdf",
            width  = 8,
            height = 7)
  # Set up plot area
  par(mfrow   = c(2, 1),
      mar     = c(0, 0, 0, 0),
      oma     = c(5, 3, 0, 0) + 0.25,
      mgp     = c(1.5, 0.5, 0),
      pty     = "m",
      lheight = 0.5)
  # plot figure
  fig_maxPD(metric)
  # end plotting
  dev.off()
}

{
  # plot multiple distance metrics
  cairo_pdf("fig/figS3-pairwise_all.pdf",
            width  = 15,
            height = 5)
  # set-up plot area
  par(mfcol = c(2, 5),
      pty = "m",
      mar = c(0, 2, 0, 1),
      oma = c(5.5, 2, 4, 0) + 0.25,
      mgp = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  fig_pdAll()
  # end plotting
  dev.off()
}

# 1.2.1. PD rarefaction #
# --------------------- #

pblapply(bin_data, function(scheme) {
  # length of rarefactions to do
  n_rar <- seq(2, length(scheme$bin_data[[1]]))

  # taxa sampled
  rar_tax_samp <- lapply(n_rar, function(n_samp) {
                        replicate(n = 500,
                                  sample(bin_data$epochs$bin_data[[1]],
                                         size = n_samp,
                                         replace = FALSE))
                   })

  mpd_rar <- lapply(rar_tax_samp, function(rar) {
    apply(rar, 2, function(repl) {
      mean(dist_data[[4]][repl,repl], na.rm = TRUE)
    }) %>%
      sort()
  })
  wmpd_rar <- lapply(rar_tax_samp, function (rar) {
    apply(rar, 2, function(repl) {
      sum(dist_data[[4]][repl, repl] * dist_data[[5]][repl, repl], na.rm = TRUE) /
        sum(dist_data[[5]][repl, repl], na.rm = TRUE)
    }) %>%
      sort()
  })

  dat <- data.frame(n = n_rar,
             mpd = sapply(mpd_rar, mean),
             mpd_lower = sapply(mpd_rar, function(x) x[length(x) * 0.025]),
             mpd_upper = sapply(mpd_rar, function(x) x[length(x) * 0.975]),
             wmpd = sapply(wmpd_rar, mean),
             wmpd_lower = sapply(wmpd_rar, function(x) x[length(x) * 0.025]),
             wmpd_upper = sapply(wmpd_rar, function(x) x[length(x) * 0.975]))



      })




# 1.3. PCo analyses #
# ----------------- #

# trim incomparable taxa from distance matrices
trimmed_data <- lapply(dist_data, Claddis::TrimMorphDistMatrix)

# remove incomparable taxa from GED matrix for consistency, and outgroup
# _Hupehsuchus nanchangenesis_
taxa_remove <- c(trimmed_data$raw.dist.matrix$removed.taxa,
                 outgroup)
trimmed_data <- lapply(trimmed_data, "[[", 1)[-5]
trimmed_data <- lapply(trimmed_data, function(mat) {
                  mat <- mat[!rownames(mat) %in% taxa_remove,
                             !colnames(mat) %in% taxa_remove]
})

# rename to shorted names
names(trimmed_data) <- plot_names$dist

# perform a principal coordinates analysis on each distance matrix
# first without Caillez (1983) correction, second without
# each list item has a `dist` and `corr` element indicating the distance matrix
# and whether Caillez correction was used for identification downstream
# also does Pearson correlation between distance and per-axes PCo data
pco_data <- c(lapply(seq_along(1:4), function (mat) {
                # get matrix
                mat <- trimmed_data[mat]

                # pco analysis
                pco <- cmdscale(mat[[1]], k = nrow(mat[[1]]) - 1, add = FALSE)
                # per-axis variance for scree plots
                scree <- scree(pco)

                # do PCo axes correlation tests
                pc_cor <- pcCorrelation(pco, mat[[1]])

                # return a list
                return(list(dist   = names(mat),
                            corr   = "uncorrected",
                            pco    = pco,
                            scree  = scree,
                            pc_cor = pc_cor))
              }),
              lapply(seq_along(1:4), function (mat) {
                # get matrix
                mat <- trimmed_data[mat]

                # pco analysis
                pco <- cmdscale(mat[[1]], k = nrow(mat[[1]]) - 1, add = TRUE)$points
                # per-axis variance for scree plots
                scree <- scree(pco)

                # do PCo axes correlation tests
                pc_cor <- pcCorrelation(pco, mat[[1]])

                # return a list
                return(list(dist  = names(mat),
                            corr  = "Caillez-corrected",
                            pco   = pco,
                            scree = scree,
                            pc_cor = pc_cor))
              }))
saveRDS(pco_data, file = "output/pco_data.rds")




# get scree values
scree_data <- lapply(pco_data, function (pco) {
                return(list(scree  = pco[["scree"]],
                            pc_cor = pco[["pc_cor"]]))
                })

{
  # scree plots of PCo data
  cairo_pdf("fig/figS4-scree_plot.pdf",
            width  = 8,
            height = 5)
  # set-up plot area
  par(mfrow = c(2, 4),
      mar   = c(0, 0, 0, 0),
      oma   = c(4, 5, 4, 4),
      mgp   = c(3, 0.5, 0))
  # plot figure
  fig_scree()
  # end plotting
  dev.off()
}


# 1.4. Morphospace #
# ---------------- #

# use uncorrected PCo from MAX distances (i.e. element 4)
pc_matrix <- 4

# PCo axes to use
PC <- list(x = 1,
           y = c(2, 3))

# get colours and points for plotting
epoch_col <- bin_colours$epochs
epoch_pch <- rep(c(15, 16, 17, 18), 2)

{
  # plot morphospace for MAX distances
  cairo_pdf("fig/fig2-max_morphospace.pdf",
            width  = 5,
            height = 7)
  # set-up plot area
  par(mfcol   = c(2, 1),
      mar     = rep(0, 4),
      oma     = c(6, 3, 0, 0) + 0.25,
      mgp     = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  fig_xyMorphospace()
  # end plotting
  dev.off()
}

{
  # plot morphospace for all PCo analyses
  cairo_pdf("fig/figS7-all_morphospace.pdf",
            width   = 12,
            height  = 5,
            onefile = TRUE)
  # set-up plot area
  par(mfcol = c(2, 4),
      mar = c(0, 3, 0, 1),
      oma = c(5, 0, 4, 0) + 0.25,
      mgp = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  fig_xyMorphospaceAll()
  # end plotting
  dev.off()
}


# 1.5. dispRity analyses #
# ---------------------- #

# export packages and variables to cluster for parallelisation
clusterCall(clus, function () {
  # required libraries
  library(dispRity)
  library(magrittr)
})
clusterExport(clus, c("bins", "subsetBootstrap", "taxon_ages", "bin_data"))

# subset by bin for each binning scheme then bootstrap and rarefy
# creates a list with index elements for binning scheme
subset_boot <- lapply(seq_along(bins), function (scheme) {
                  # apply binning into each `pco_data` list element
                  lapply(pco_data, function (pco) {
                    # combine into list
                    return(c(pco,
                             bins = names(bins[scheme])))
                  })
                }) %>%
                unlist(recursive = FALSE) %>%
                pblapply(cl = clus, function (run) {
                  # create binning vector for dispRity function
                  # bin_splits <- c(bins[[run$bins]]$max_age,
                  #                 min(bins[[run$bins]]$min_age))
                  bin_splits <- bin_data[[run$bins]]$bin_data
                  # remove outgroup
                  bin_splits <- bin_splits[names(bin_splits) != "outgroup"]

                  # keep only taxa in the matrix
                  bin_splits <- lapply(bin_splits, function (bin) {
                                  bin[bin %in% rownames(run$pco)]
                                })

                  # subset & bootstrap data
                  sub <- subsetBootstrap(run$pco,
                                         bin_splits)

                  # return data:
                  # combine list from pco_data with bin scheme and
                  # dispRity object
                  return(c(run,
                           list(subboot = sub)))
                })

# do disparity analyses; parallelised for speed
# assign disparity metric then do calculations; allows for per-run parallelisation
disp_calc <- lapply(seq_along(plot_names$disp_metrics), function (metric) {
                # apply metrics into each `subset_boot` list element
                lapply(subset_boot, function (run) {
                  # combine into list
                  return(c(run,
                           metric = names(plot_names$disp_metrics[metric])))
                })
              }) %>%
              unlist(recursive = FALSE) %>%
              pblapply(cl = clus, function (run) {
                # disparity calculation
                disp <- dispRity(run$subboot,
                                 metric = plot_names$disp_metrics[[run$metric]][["func"]])

                # return data:
                # combine list from subset_boot with disparity object
                return(c(run,
                         list(disp_calc = disp)))
              })
saveRDS(disp_calc, "output/disp_calc.rds")
# send notification email when finished
sendEmailNotification(subject = "Notification: dispRity success",
                      body    = "dispRity analyses completed successfully. :)")

# summarise dispRity object to get mean disparity and 95% confidence intervals
disparity <- pblapply(disp_calc, function (run) {
                # summarise dispRity object
                disp <- summary(run$disp_calc,
                                quantiles = 95,
                                cent.tend = mean,
                                digits = 5)

                # get other metrics for t-tests
                summ_test <- lapply(seq_along(run$disp_calc$disparity), function (nbin) {
                               # number of run
                               bin <- run$disp_calc$disparity[[nbin]]

                               # data.frame
                               data.frame(bs.mean = mean(bin[[2]]),
                                          st_d    = sd(bin[[2]]),
                                          samp    = length(run$disp_calc$subsets[[nbin]]$elements),
                                          `2.5%`  = quantile(bin[[2]],
                                                            probs = 0.025,
                                                            na.rm = TRUE),
                                          `97.5%` = quantile(bin[[2]],
                                                             probs = 0.975,
                                                             na.rm = TRUE),
                                          check.names = FALSE)
                             }) %>%
                             do.call(rbind, .)

                # change the row names to epoch names (for epochs only)
                if (run$bins == "epochs") {
                  rownames(summ_test) <- rownames(bin_data$epochs$bins)
                } else {
                  rownames(summ_test) <- rownames(bin_data$tenMa_bJ$bins)
                }

                # return data:
                # combine list from disp_calc with summarised disparity
                return(list(dist      = run$dist,
                            corr      = run$corr,
                            bins      = run$bins,
                            metric    = run$metric,
                            disp_summ = disp,
                            summ_test = summ_test))
             })
saveRDS(disparity, file = "output/disparity.rds")

{
  # plot disparity curves
  cairo_pdf("fig/figS5-disparity_time.pdf",
            width   = 12,
            height  = 5,
            onefile = TRUE)
  # set-up plot area
  # TODO: check scaling etc
  par(mfcol   = c(2, 4),
      mar     = c(0, 2, 0, 2),
      oma     = c(5.5, 2, 4, 2),
      mgp     = c(1.5, 0.5, 0),
      lheight = 0.5)
  # plot figure
  fig_disparityAll()
  # stop plotting
  dev.off()
}


# 1.6. Rarefaction curves #
# ----------------------- #

{
  # plot rarefaction curves
  cairo_pdf("fig/figS6-rarefaction_curves.pdf",
            width = 10,
            height = 10,
            onefile = TRUE)
  # set-up plot area
  par(oma = c(3, 4, 2, 1),
      mar = c(1, 1, 3, 1),
      mgp = c(1.5, 0.5, 0))
  # plot figure
  fig_rarefactionCurves()
  # stop plotting
  dev.off()
}

# stop cluster
stopCluster(clus)
