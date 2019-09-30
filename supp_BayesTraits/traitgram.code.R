#1. Preparation ####
rm(list=ls())

# wd <- "C:/Users/tstub/Dropbox/ichthy-macroevo/BayesTraits_results"

wd <- getwd() 
wd

#Check which operating system we are using and set the path accordingly
#We only check for Windows and default otherwise to Linux/Unix
#Before setting the path we check whether we are running in
#interactive mode or not - if not, the path to the invoked script
#is assumed to represent the working directory (works only on Linux)
if(interactive() == TRUE){
  if(.Platform$OS.type == "windows"){
    myOS <- "Windows"
    } else{
    myOS <- "NOT_Windows"
  }
} else{
  myOS <- "NOT_Windows"
  setwd(system("pwd", intern = T)) #set path in non-interactive mode (i.e. when using Rscript) - works only on Linux
}

# Which clade should we focus on
myclade = "Ichthyosaurs" #Switch between different clades

# 2. Install and load packages####

library(ggplot2)
library(ape)
library(geiger)
library(strap)
library(paleotree)
library(phytools)
library(nlme)
library(caper)
library(coda)
library(rio)
library(lattice) #necessary for some of the convergence plots
library(foreach) #load foreach library (needed to parallelize some of our tasks)
library(parallel) #load parallel library (needed to parallelize some of our tasks)
library(doParallel) #load doParallel library (needed to parallelize some of our tasks)

sessionInfo()

#############################################################

#2. Import data####
Ichthyosaur.tree <- read.tree("sample_trees.tre")
Ichthyosaur.tree

# plot(ladderize(Ichthyosaur.tree[[1]]),cex=0.5)
# plot(ladderize(Ichthyosaur.tree[[100]]),cex=0.5)

#Import the .txt file used to time-scale the phylogeny
Ichthyosaur.ages <-read.csv("ichthyosaur_occurrences.csv", header=T, row.names=1)
Ichthyosaur.ages

####3. Time-scale the phylogeny####
####Preparations####

name.check(Ichthyosaur.tree[[1]],Ichthyosaur.ages) # Check if all species names are present in both
myclade.ages <- Ichthyosaur.ages[! rownames(Ichthyosaur.ages) %in% name.check(Ichthyosaur.tree[[1]], Ichthyosaur.ages)$data_not_tree,]
name.check(Ichthyosaur.tree[[1]],myclade.ages) #Check if all species names are present in both

###### Need to resolve base of phylogeny here
Ichthyosaur.tree.rooted <- lapply(Ichthyosaur.tree, root, outgroup = "Hupehsuchus_nanchangensis", resolve.root = TRUE)
class(Ichthyosaur.tree.rooted) <- "multiPhylo"
Ichthyosaur.tree.rooted
Ichthyosaur.tree.rooted[[1]]

Ichthyosaur.tree.rooted.no.edges <-lapply(Ichthyosaur.tree.rooted,function(x){x$edge.length<-NULL; x})
class(Ichthyosaur.tree.rooted.no.edges)<-"multiPhylo"
Ichthyosaur.tree.rooted.no.edges
Ichthyosaur.tree.rooted.no.edges[[1]]

# plot(ladderize(Ichthyosaur.tree.rooted.no.edges[[1]]),cex=0.5)
# plot(ladderize(Ichthyosaur.tree.rooted.no.edges[[1]]),cex=0.5)
# comp.trees <- cophylo(Ichthyosaur.tree.rooted.no.edges[[1]],Ichthyosaur.tree.rooted.no.edges[[50]])
# plot(comp.trees)

########################################################################
########################################################################
########################################################################
########################################################################

raw.mid_age <- data.frame(matrix(nrow=length(myclade.ages[,1]),ncol=2))
raw.mid_age[,1] <- as.matrix(myclade.ages[, "LAD"])
raw.mid_age[,2] <- as.matrix(myclade.ages[, "LAD"])
row.names(raw.mid_age) <- row.names(myclade.ages)
raw.mid_age
dim(raw.mid_age)

# DATE THE PHYLOGENY USING THE NODE DATA AND EQUAL METHOD FOR ALL OTHER NODES
time.tree <- timePaleoPhy(Ichthyosaur.tree.rooted.no.edges[[1]], raw.mid_age, type = "equal", vartime=10, ntrees=1, randres=FALSE, dateTreatment = "firstLast", node.mins = NULL)

plot(ladderize(time.tree), cex=0.4, no.margin = TRUE)

#4  Import the body size data ## 
df <- read.csv("ichthyosaur_data_BAYESTRAITS.csv", row.names=1, header=T)
 
bodySize <- matrix(df$Skulllength.mm) # chose data of interest
rownames(bodySize) <- rownames(df)
nms.id <- complete.cases(bodySize)
bodySize.pruned <- as.matrix(bodySize[nms.id,]) #Remove all the taxa for which we have no data from our matrix
bodySize.pruned
dim(bodySize.pruned)

####### Log-transform our data ########
bodySize.pruned <- log10(bodySize.pruned)
bodySize.pruned

droptips_trees <- name.check(time.tree, bodySize.pruned) #Find taxa to prune using one random tree
droptips_trees 
bodySize.pruned <- bodySize.pruned[!rownames(bodySize.pruned) %in% droptips_trees$data_not_tree,]
bodySize.pruned
length(bodySize.pruned)

## Now drop the tips from ALL our timescaled and randomly resolved trees
time.tree <- drop.tip(time.tree,droptips_trees$tree_not_data)
time.tree
plot(ladderize(time.tree), cex=0.4)
time.tree$root.time

###############################################################################################
################################################################################################
################################################################################################

plot_tree <- ladderize(time.tree) # select tree
plot(plot_tree)

plot_data <-bodySize.pruned  # select data

trait <- plot_data
traitname <- "Skull length (log10 mm)"

# rates.results <- read.table("edge.rates.txt")
# palette <- colorRampPalette(colors=c("darkblue","orange", "red"))
# cols <- palette(98); names(cols)<-1:98

treeheight<-max(vcv(plot_tree))
y <- c(1.7,(max(plot_data)))
# y <- c(round(min(crop.data[, plot_metric2])-(0.1*max(crop.data[, plot_metric2]))):round(max(crop.data[, plot_metric2])+(0.1*max(crop.data[, plot_metric2]))))
x <-seq(0,treeheight,treeheight)
line1<-plot_tree$root.time-250
line2<-plot_tree$root.time-200
line3<-plot_tree$root.time-150
line4<-plot_tree$root.time-100
xat<-c(line1, line2, line3, line4)
xlabels<-c("250","200","150", "100")

par(mai=c(0.5,0.5,0.2,0.2),mgp=c(1.5,0.5,0), las=1)
plot(x,y, type="n", axes=FALSE,ylab=traitname, xlab="Mya", cex.lab=1)
abline(v=c(pbound, jbound,cbound),col="grey",lwd=1.5,lty=2)
axis(2, cex.axis=0.7)
axis(side = 1, at = xat, labels = xlabels, tck=-.02, cex.lab=0.7,cex.axis=0.7)
box()
try <- phenogram(plot_tree, trait,fsize=0,spread.labels=F, colors="gray50", add=TRUE)

# save plot
dev.copy(pdf,'skull.size.traitgram.pdf', width = 7, height = 5)
dev.off()

######
plot(x,y, type="n", axes=FALSE,ylab=traitname, xlab="Mya", cex.lab=1)
phenogram(plot_tree, trait,fsize=1,spread.labels=T, colors="dodgerblue2")
# save plot
dev.copy(pdf,'traitgram.label.pdf', width = 15, height = 10)
dev.off()

colors=colors.branches

polygon (c(255, 245, 245, 255), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(235, 225, 225, 235), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(215, 205, 205, 215), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(195, 185, 185, 195), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(175, 165, 165, 175), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(155, 145, 145, 155), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(135, 125, 125, 135), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))
polygon (c(115, 105, 105, 115), y = c(-1000, -1000, 3000, 3000), border = NA, col = rgb (0, 0, 0, 0.03))

box()
mtext("Relative evolutionary rates", 2, line = 3, cex = 0.8, font = 2)
mtext("Time (Ma)", 1, line = 3.5, cex = 0.8, font = 2)
axis(side=2)


