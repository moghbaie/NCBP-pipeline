library(data.table)
library(stringr)
library(R.utils)
require(methods)

setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd('..')

source("src/MaxQuant.R")
source("src/Preprocessing/preprocessing.R")

source("src/Imputation/impute.R")

source('src/Stat.analysis/ANOVA.R')
source('src/Stat.analysis/select_significant.R')
source('src/Stat.analysis/Plot_volcano_plots.R')
source('src/Stat.analysis/CompareSeeds.R')
source('src/Stat.analysis/ANOVA_vote.R')
source('src/Stat.analysis/Plot_volcano_plots.R')
source("src/Stat.analysis/control_FP_FN.R")


source("src/Clusterization/distances.R")
source("src/Clusterization/make_clusters.R")

source("src/QC/QC1.R")

source("src/Visualization/protein_copy_number.R")
source("src/Visualization/plot_disco.R")
source("src/Visualization/plot_wave.R")
source("src/Visualization/plot_heatmap.R")
source("src/Visualization/plot_behaviour.R")
source("src/Visualization/MDS.R")


# PARAMETERS
recalculate.seeds = F
N.seeds = 100

# QC1
PerformQC()


for (pattern in c("LFQ intensity", "iBAQ")) {
  max.quant <- BuildMaxQuant("./data/txt_noGA/proteinGroups.txt",
                             pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList)
  max.quant$MakePrettyExperimentNames()
  
  # initial data preparation
  max.quant$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky() 
  
  if (pattern == "LFQ intensity" & recalculate.seeds) {
    RunNtimes(N.seeds, max.quant, n.threads = 7)
    PlotNseedsHeatmap(N.seeds, max.quant)
    PlotAmountOfSignificantHeatmap(N.seeds, max.quant)
  }
  # imputation
  max.quant$ImputeAll(id.column = kIdColumnName, sep.method = "mnar.1.nonzero",
                      method.mar = "mar.stat.6", 
                      method.mnar = 'mnar.stat.1', sigma.const = 4,
                      save.result.table = T, seed = 46)
  # normalization
  max.quant$
    PowTwo()
  raw.intensities.range <- FindIntensityRange(max.quant)
  max.quant$
    Normalize()$
    Rescale(range = raw.intensities.range, by.coef = T)$
    LogTricky()
  
  if (pattern == "LFQ intensity") {
    max.quant$AnovaVote(sign.cutoff = 60, N.seeds)
    
    PlotVolcanoPlots(max.quant, mode = 'outside')
    # Control FP & FN
    max.quant$ControlFPFN(passed.anova = F, save.fp.0.1 = T, mode = "vote")
    max.quant$ControlFPFN(passed.anova = T, save.fp.0.1 = T, mode = "vote")
    
    PlotDiscoPlots(max.quant, mode = 'target')
    PlotDiscoPlots(max.quant, mode = 'buffer')
    # select significant proteins
    max.quant.sign <- max.quant$MakeObjectOnlySignificant(mode = 'vote')
    
    PlotWave(max.quant.sign, n.anovas.passed = c(1, 2)) 
  } else if (pattern == 'iBAQ') {
    max.quant.sign <- max.quant$MakeObjectOnlySignificant(mode = 'vote')
  }
  
  max.quant.sign$
    CountDistanceMatrix(metrics = "euclidean", all.targets = T)$
    MakeClusters(n.of.clusters = 5)
  
  PlotHeatMap(max.quant.sign, color = 'colorful', buffers.mode = "reorder", reorder.by.dist = T,
              make.table.of.distances = F, control.mode = 'subtract',  metrics = 'euclidean', 
              recalculate.dist = F)
  PlotHeatMap(max.quant.sign, color = 'colorful', buffers.mode = "only shared", reorder.by.dist = T,
              make.table.of.distances = F, control.mode = 'subtract',  metrics = 'euclidean', 
              recalculate.dist = F)
  PlotDendrosWithBehaviour(max.quant.sign, n.of.clusters = 5) 
  PlotMDSWithArrows(max.quant.sign) 
}

