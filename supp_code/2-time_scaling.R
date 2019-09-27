#!/usr/bin/env Rscript

# 2. Time-scaling trees #
# ===================== #

# Script from Moon & Stubbs <Title>, <doi>

# required libraries
library(Claddis)
library(doParallel)
library(magrittr)
library(parallel)
library(pbapply)
library(strap)

# required scripts
source("functions/functions_7.R")
source("functions/scaling_functions.R")
source("functions/sendEmailNotification.R") # sends email; requires set-up

# for parallelisation
no_cores <- detectCores() - 1
clus <- makeCluster(no_cores, outfile = "")
registerDoParallel(clus)


# 2.1. Gather data #
# ---------------- #

# load tree, age, and cladistic data
trees <- ape::read.tree("data/sample_trees.tre")[1:120]
nexus_data <- ReadMorphNexus("data/matrix.nex")
full_ages <- read.table("data/ichthyosaur_occurrences.tsv",
                        sep       = "\t",
                        header    = TRUE,
                        row.names = 1)
ages <- full_ages[match(trees[[1]]$tip.label, rownames(full_ages)), ]

# root, resolve, and ladderize trees
root_trees <- lapply(trees, function(x) {
                ape::root(x,
                          outgroup = "Hupehsuchus_nanchangensis",
                          resolve.root = TRUE) %>%
                ladderize
              })
class(root_trees) <- "multiPhylo"

# write `root_trees` to Newick file and re-read:
# the Hedman method sometimes produces negative-length branches, but this seems
# to resolve that
write.tree(root_trees, file = "output/root_trees.tre")
root_trees <- read.tree("output/root_trees.tre")

# load outgroup ages 
outgroup_ages <- read.table("data/outgroup_occurrences.tsv",
                            sep       = "\t",
                            header    = TRUE,
                            row.names = 1)


# 2.2. Minimum branch length time-scaling #
# --------------------------------------- #

# export paleotree package to cluster
clusterCall(clus, function () {
  library(paleotree)
})
clusterExport(clus, c("ages"))

# time-scale trees using a MBL of 1 Ma, and minMax observations; 10 trees from
# each run are returned
mbl_ttrees <- pblapply(root_trees, cl = clus, function (tree) {
                # do MBL tree scaling
                tt <- timePaleoPhy(tree,
                                   ages,
                                   type    = "mbl",
                                   vartime = 1,
                                   ntrees  = 10,
                                   dateTreatment = "minMax",
                                   noisyDrop = FALSE)

                # add MBL ID element to TS-tree
                tt <- lapply(tt, function (ttree) {
                        return(list(tree    = ttree,
                                    scaling = "mbl"))
                      })

                # return time-scaled trees
                return(tt)
              }) %>%
              unlist(recursive = FALSE)


# 2.3. Hedman (2010) time-scaling #
# ------------------------------- #

# time_scale trees using Hedman (2010) method; incorporates functions by Lloyd
# from script `functions_7.R`

# pull occurrence dates for each tree from a uniform distribution
time_uniform <- lapply(root_trees, function (tree) {
                  # occurrences for the ingroup
                  ingroup <- timeHedman(ages)

                  # outgroup occurrences
                  outgroup <- timeHedman(outgroup_ages)

                  # return list
                  return(list(ingroup  = ingroup,
                              outgroup = outgroup))
})

# Hedman time-scaling using outgroup mid-ages; failed runs are removed
{
  hedman_ttrees <- list()
  hedman_ttrees <- foreach (tree = seq_along(root_trees), .errorhandling = "remove") %dopar% {
    # Hedman scaling
    out <- Hedman.tree.dates(root_trees[[tree]],
                             time_uniform[[tree]]$ingroup,
                             outgroup_ages$MID,
                             t0           = 400,
                             resolution   = 10000,
                             conservative = TRUE)

    # return tree and scaling method
    hedman_ttrees[[tree]] <- list(tree    = out$tree,
                                  scaling = "Hedman")
  }
  # send notification email when finished
  sendEmailNotification(subject = "Notification: Hedman success",
                        body    = "Hedman tree-scaling completed successfully. :)")
}

# save trees and write to files
saveRDS(hedman_ttrees, file = "output/hedman_ttrees.rds")

# write trees to Newick and Nexus files
hedman_write <- lapply(hedman_ttrees, "[[", "tree")
class(hedman_write) <- "multiPhylo"
write.tree(hedman_write,
           file = "data/hedman_ttrees.tre")
write.nexus(hedman_write,
            file = "data/hedman_ttrees.nex",
            translate = TRUE)

# combine all trees
time_trees <- c(mbl_ttrees, hedman_ttrees)

# save time_trees object
saveRDS(time_trees, file = "output/time_trees.rds")


# stop cluster
stopCluster(clus)