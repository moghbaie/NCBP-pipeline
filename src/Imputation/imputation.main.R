library(data.table)
library(stringr)
require(pbmcapply)
require(ggplot2)
require(plyr)

source("src/MaxQuant.R")
source("src/Preprocessing/preprocessing.R")

source("lib/constants.R")
source("lib/helper_functions.R")

source("src/Imputation/impute.R")
source('src/Imputation/replicates.correlation.R')
source("src/Imputation/set.null.and.impute.R")
source("src/Imputation/rmse.R")
source("src/Imputation/auc.R")
source("src/Imputation/correlation_matrix.R")

source("src/Stat.analysis/control_FP_FN.R")
source("src/Stat.analysis/ANOVA.R")
source("src/Stat.analysis/select_significant.R")


patterns <- c("LFQ intensity", "iBAQ")
for (pattern in patterns) {
  CheckExperimentsCorr(path = "./data/txt_noGA/proteinGroups.txt", pattern = pattern,
                       sep.method = "mnar.1.nonzero", method.mar = "mar.stat.6", 
                       method.mnar = "mnar.stat.1")
  CompareWithPerseus(path.to.perseus.output = "./data/NCBP3(Columbia) Perseus output")
  CheckReplicatesCorrelationAfterImputation()
  
  max.quant <- BuildMaxQuant("./data/txt_noGA/proteinGroups.txt", pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList)
  max.quant$MakePrettyExperimentNames()
  max.quant$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky() 
  
  CountRMSES(max.quant, n.threads = 2, n.iterations = 2)
  PlotROCcurves(max.quant)
  MakeAUCsConfIntervals(max.quant, n.threads = 2, n.iterations = 2)  
}



