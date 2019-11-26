list.of.packages <-
  c(
    "ggplot2",
    "data.table",
    "stringr",
    "R.utils",
    "R6",
    "scales",
    "dplyr",
    "digest",
    "dendextend",
    "plotly",
    "ape",
    "vegan",
    "proxy",
    "reshape2",
    "bit64",
    "ggthemes",
    "reshape",
    "PTXQC",
    "ade4",
    "pbmcapply",
    "clValid",
    "gplots", 
    "cba",
    "ggrepel",
    "circlize",
    "usedist",
    "foreach",
    "ggnewscale",
    'doParallel'
  )
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("UniProt.ws")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
