library(R.utils)
source("src/QC/proteins.presence.R")
source("src/QC/contaminants.R")
source("src/QC/plot.intervals.of.intensities.R")
source("src/QC/compare.dist.mtx.R")
source("src/QC/explore.close.isoforms.R")
source("src/QC/PTXQC.R")




PerformQC <- function(){
  QCReport('data/txt_noGA/')

  ## distribution of intensities (to exclude data)
  pattern <- 'iBAQ'
  max.quant <- BuildMaxQuant("./data/txt_LFQ_no_match_groups_NCBP3x2_new/proteinGroups.txt", pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList, new.data = F)
  max.quant$
    MakePrettyExperimentNames()$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky()
  PlotTargetIntensityHist(max.quant, kTargetProteins, new.data = F)

  ## TIC vs contaminants (to exclude data)
  PlotTICvsContaminants(path.to.mq = "data/txt_LFQ_no_match_groups_NCBP3x2_new/")
  
  pattern = "LFQ intensity"
  max.quant <- BuildMaxQuant("./data/txt_LFQ_no_match_groups_NCBP3x2_new/proteinGroups.txt", pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList, new.data = F)
  max.quant$
    MakePrettyExperimentNames()$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky()
  max.quant$ImputeAll(id.column = kIdColumnName, sep.method = "mnar.1.nonzero",
              method.mar = "mar.stat.6", 
              method.mnar = 'mnar.stat.1', sigma.const = 4,
              save.result.table = F, seed = 46)
  
  max.quant$
    PowTwo()
  raw.intensities.range <- FindIntensityRange(max.quant)
  max.quant$
    Normalize()$
    Rescale(range = raw.intensities.range, by.coef = T)$
    LogTricky()
  
  ## Mantel test (to exclude data)
  max.quant$PerformANOVA(do.plot = F, save.result.table = F)
  max.quant.sign <- max.quant$MakeObjectOnlySignificant(mode = 'single')
  max.quant.sign$
    CountDistanceMatrix(metrics = "euclidean", all.targets = T)
  distance.mtx1 <- max.quant.sign$GetDistanceMatrix()
  targets <- c("Rockefeller_NCBP3")
  max.quant$FilterColumns(pattern = paste(targets, collapse = '|'))
  max.quant.sign <- max.quant$MakeObjectOnlySignificant(mode = 'single')
  max.quant.sign$
    CountDistanceMatrix(metrics = "euclidean", all.targets = T)
  distance.mtx2 <- max.quant.sign$GetDistanceMatrix()
  CompareDistanceMtxs(distance.mtx1, distance.mtx2)

  ## explore MAGOH/B
  ExploreCloseIsoforms(path.to.mq = "data/txt_noGA/", isoform1 = "P61326", isoform2 = "Q96A72", old.data = F)
  
  
  
  ## Normalization check
  pattern = "LFQ intensity"
  max.quant <- BuildMaxQuant("./data/txt_noGA/proteinGroups.txt", pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList, new.data = T)
  max.quant$
    MakePrettyExperimentNames()$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky()$
    ImputeAll(id.column = kIdColumnName, sep.method = "mnar.1.nonzero",
              method.mar = "mar.stat.6", 
              method.mnar = 'mnar.stat.1', sigma.const = 4,
              save.result.table = F, seed = 1)
  PlotIntervalOfIntensities(max.quant, proteins.to.plot = c(kTargetProteins, kGFPprotein))
  max.quant$
    PowTwo()
  raw.intensities.range <- FindIntensityRange(max.quant)
  max.quant$
    Normalize()$
    Rescale(range = raw.intensities.range, by.coef = T)$
    LogTricky()
  PlotIntervalOfIntensities(max.quant, proteins.to.plot = c(kTargetProteins, kGFPprotein))
}

