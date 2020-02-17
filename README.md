[![DOI](https://zenodo.org/badge/211271780.svg)](https://zenodo.org/badge/latestdoi/211271780)

---
title: Early high rates and disparity in the evolution of ichthyosaurs 
author:
- name: Benjamin C. Moon
  affiliation: University of Bristol
  email: benjamin.moon@bristol.ac.uk
- name: Thomas L. Stubbs
  affiliation: University of Bristol
  email:  tom.stubbs@bristol.ac.uk
date: 2019-12-09
journal: Communications Biology
url: https://www.nature.com/articles/s42003-020-0779-6
doi: https://doi.org/10.1038/s42003-020-0779-6
date_published: 2020-02-13
...

# ichthyosaur-macroevolution #

Data and code used in our study ‘Early high rates and disparity in the
evolution of ichthyosaurs’, published as

> Moon BC & Stubbs TL Early high rates and disparity in the evolution of
ichthyosaurs. Communications Biology, 3, 68 doi:10.1038/s42003-020-0779-6

# Contents #

1. supplement-ichthyosaur_macroevolution: supplementary document PDF file with
   supplementary methods, results, figures, tables, and references.
2. supp_figures: folder of supplementary figures included in the supplementary
   document.
3. supp_tables: folder of supplementary  tables included in the supplementary
   document.
4. supp_code: supplementary R code to run most analyses.
5. supp_BayesTraits: supplementary code to run the continuous rates analyses in
   BayesTraits.

## Repository layout ##

    ichthyosaur-macroevolution
    ├ LICENSE
    ├ README.md
    ├ supplement-ichthyosaur_macroevolution.tex
    ├ supplement-ichthyosaur_macroevolution.pdf
    ├ supplement-ichthyosaur_macroevolution.bib
    ├ version.dat
    ├ supp_figures/
    ├ supp_tables/
    ├ supp_code/
    │  ├ 0-Moon_Stubbs-run_analyses.R
    │  ├ 1A-disparity_analyses.R
    │  ├ 1B-disparity_stats.R
    │  ├ 2-time_scaling.R
    │  ├ 3-rates_analyses.R
    │  ├ data/
    │  ├ fig/
    │  ├ functions/
    │  ├ output/
    │  └ tbl/
    └ supp_BayesTraits/
       ├ ichthyosaur_data_BAYESTRAITS.csv
       ├ Moon_Plot_Functions.R
       ├ traitgram.code.R
       ├ BayesTraits/
       └ PostProc/

# Running the code #

## Compiling the supplement ##

This uses recent versions of `lualatex` (compiled with version 1.10.0) and
`latexmk` (version 4.65). Run the following line to compile

    > latexmk supplement-ichthyosaur_macroevolution

## R code analyses ##

Most of the analyses are run in the R statistical environment; these were
tested using R version 3.6.0. Scripts are contained in the folder `supp_code/`
All analyses will proceed continuously from the first script from the shell

    > Rscript 0-Moon_Stubbs-run_analyses.R

or within R

    r$> source("0-Moon_Stubbs-run_analyses.R")

This includes installation of the required packages at the beginning of the
first script. Refactoring of Claddis took place during the analyses, so we use
a previous version (commit 4f7f1bf); the original version of this code will not
function with a later commit. Email notifications are sent using the packages
`mailR` and password storage with `keyringr`. Vignettes for the set-up of these
can be found at <https://github.com/rpremraj/mailR> and
<https://bit.ly/2AHUj0A> (accessed 2019-09-27).

## BayesTraits analyses ##

Individual command scripts and associated trees for 100 runs are included as
`.cmd` and `.nex` files respectively in the subfolder
`supp_BayesTraits/BayesTraits`. Run these in BayesTraits (version 2.0.2) from
<http://www.evolution.rdg.ac.uk/BayesTraitsV2.0.2.html> for automated analyses

    > BayesTraits BayesTraits/ichthyosaur.pca.scores.pruned.singletree.1.nex
    ichthyosaur_data_BAYESTRAITS.csv
    < BayesTraits/BayesTraits_hetRates_Script.1.cmd

Homogeneous and heterogeneous scripts are separated as `homRates` or
`hetRates`. Results tables are included in the folder
`supp_BayesTraits/PostProc`. Additional plot are created using the scripts
`traitgram.code.R` and `Moon_plot_Functions.R`.
