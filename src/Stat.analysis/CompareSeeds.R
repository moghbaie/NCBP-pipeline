library(data.table)
library(stringr)
library(R.utils)
library(reshape)
library(reshape2)
library(plyr)
library(viridis)
library(parallel)
library(pbmcapply)
library(doSNOW)

RunNtimes <- function(N.seeds, max.quant.obj, n.threads){
  dir.create(file.path('out/QC/FPFN_seeds/'), showWarnings = FALSE)
  dir.create(file.path('out/QC/pvals_logfc_seeds/'), showWarnings = FALSE)
  cl <- makeCluster(n.threads, type = "SOCK")
  clusterExport(cl, c("ImputeCaseAndControl", "Impute", 
                      "SeparateImputeMARMNARMerge", "SeparateImputeByNNonZeroMerge",
                      "SeparateByNNonZero", "ImputeMNAR", "ImputeMAR", "MergeAll", 
                      "kIdColumnName", "FindIntensityRange", "CountMeanSDNoOutliers",
                      "SeparateImputeOneColForMNAR", "MakeDate", "kGroupOrder", "kNoGroupName"))
  registerDoSNOW(cl)
  max.quant.tmp <- max.quant.obj$clone()
  foreach (seed. = c(1:N.seeds), .verbose = T, .packages = c("plyr", "stringr")) %dopar% {
    library(data.table)
    max.quant.tmp$ImputeAll(id.column = 'Protein IDs', sep.method = "mnar.1.nonzero",
                            method.mar = "mar.stat.6",
                            method.mnar = 'mnar.stat.1', sigma.const = 4,
                            save.result.table = F, seed = seed.)
    max.quant.tmp$
      PowTwo()
    raw.intensities.range <- FindIntensityRange(max.quant.tmp)
    max.quant.tmp$
      Normalize(norm.method = "by_protein")$
      Rescale(range = raw.intensities.range, by.coef = T)
    max.quant.tmp$LogTricky()
    max.quant.tmp$PerformANOVA(do.plot = F, save.result.table = T, prefix = paste0('/QC/pvals_logfc_seeds/_', seed., '_'))
    max.quant.tmp$ControlFPFN(passed.anova = F, save.table = F, mode = 'single')
    saveRDS(max.quant.tmp$GetHypothesesStatus(), paste0("out/QC/FPFN_seeds/_", seed., '_', MakeDate(), "_FPFN.rds"))
  }
  stopCluster(cl)
}



PlotNseedsHeatmap <- function(N.seeds, max.quant.obj){
  metadata.experiment <- max.quant.obj$GetMetadataExperiments()
  metadata.experiment$group <- paste(metadata.experiment$location, metadata.experiment$target, metadata.experiment$buffer, sep = '_')
  metadata.protein <- max.quant.obj$GetMetadataProteins()
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Gene.names
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Protein.IDs
  metadata.protein$uniprot.gene.best <- str_split_fixed(metadata.protein$uniprot.gene.best, ';', 2)[, 1]
  results.list <- list()
  files <- list.files('./out/QC/pvals_logfc_seeds/', full.names = T)
  for(file in files){
    seed <- str_split_fixed(file, '_', 5)[, 4]
    dt <- readRDS(file)
    list.of.significant.per.column <- lapply(unique(metadata.experiment$group), function(group) {
      (dt[, paste0("log2.fold.change.", group), with = F] > 1) &
        (dt[, paste0("pvalue.adj.", group), with = F] < 0.05)
    })
    sign.m <- as.data.table(do.call("cbind", list.of.significant.per.column))
    colnames(sign.m) <- unique(metadata.experiment$group)
    sign.m$Protein.IDs <- dt$Protein.IDs
    sign.m <- replace(sign.m, is.na(sign.m), F)
    sign.m <- melt(sign.m, id.vars = 'Protein.IDs', variable.name = 'experiment', value.name = seed)
    results.list[[seed]] <- sign.m
  }
  results.dt <- merge_recurse(results.list)
  results.dt <- results.dt[rowSums(results.dt[, as.character(c(1:N.seeds)), with = F] == TRUE) > 1, ]
  results.dt <- melt(results.dt, id.vars = c('Protein.IDs', 'experiment'), variable.name = 'seed', value.name = 'is.significant')
  results.dt <- merge(results.dt, metadata.protein[, c('uniprot.gene.best', 'Protein.IDs'), with = F], by = 'Protein.IDs', all.x = T)
  
  results.list <- list()
  files <- list.files('./out/QC/FPFN_seeds/', full.names = T)
  for(file in files){
    seed <- str_split_fixed(file, '_', 4)[, 3]
    dt <- readRDS(file)
    dt <- melt(dt, id.vars = 'Protein.ID', variable.name = 'experiment', value.name = seed)
    results.list[[seed]] <- dt
  }
  FPFN.dt <- merge_recurse(results.list)
  FPFN.dt <- melt(FPFN.dt, id.vars = c('Protein.ID', 'experiment'), variable.name = 'seed', value.name = 'status')
  colnames(FPFN.dt)[colnames(FPFN.dt) == 'Protein.ID'] <- 'Protein.IDs' 
  results.dt <- merge(results.dt, FPFN.dt, by = c('Protein.IDs', 'experiment', 'seed'), all.x = T)
  target.names <- unique(metadata.protein[metadata.protein$Protein.IDs %in% grep(paste(kTargetProteins, collapse = '|'), 
                                                                                 metadata.protein$Protein.IDs, value = T) , ]$uniprot.gene.best)
  colours <- c('#F8766D', '#00BFC4')
  names(colours) <- c(FALSE, TRUE)
  
  # coloursFalses <- c('blue', 'red')
  # names(coloursFalses) <- c('FN', 'FP')
  results.dt$status.bool <- FALSE
  results.dt[results.dt$status %in% c('TP', 'TN'), ]$status.bool <- TRUE
  
  pdf(paste0('./out/QC/', MakeDate(), '_significant_proteins_heatmap.pdf'), width = 20, height = 17)
  for(group in unique(metadata.experiment$group)){
    results.group <- results.dt[results.dt$experiment == group, ]
    results.tmp <- results.group[, c('uniprot.gene.best', 'seed', 'is.significant'), with = F]
    results.tmp <- dcast(results.tmp, uniprot.gene.best~seed)
    rownames(results.tmp) <- results.tmp$uniprot.gene.best
    results.tmp <- results.tmp[, names(results.tmp) %nin% c('uniprot.gene.best')]
    ord <- sort(rowSums(results.tmp))
    protein.ord <- names(ord)
    proteins.order.dt <- as.data.frame(ord)
    proteins.order.dt$uniprot.gene.best <- rownames(proteins.order.dt)
    proteins.order.dt$ord <- proteins.order.dt$ord / N.seeds * 100
    protein.ord <- protein.ord[protein.ord %nin% target.names]
    protein.ord <- c(protein.ord, target.names)
    results.group <- merge(results.group, proteins.order.dt, by = 'uniprot.gene.best')
    results.group$Gene <- factor(results.group$uniprot.gene.best, levels = protein.ord)
    results.group$Gene2 <- paste(results.group$Gene, '_', results.group$ord, '%', sep = '')
    results.group$Gene <- mapvalues(results.group$Gene, from = unique(results.group$Gene), to = unique(results.group$Gene2))
    seed.ord <- names(sort(colSums(results.tmp)))
    results.group$seed <- factor(results.group$seed, levels = seed.ord)
    
    g <- ggplot(results.group, aes(x = seed, y = Gene, fill = is.significant)) + 
      geom_tile() +
      ggtitle(label = group) +
      scale_fill_manual(values = colours) +
      geom_point(aes(size=ifelse(status.bool, "no_dot", "dot"))) +
      scale_size_manual(values=c(dot = 1, no_dot = NA), guide = "none")
    print(g)
  }
  dev.off()
}

ReorderProteins <- function(results.dt, metadata.protein, proteins.to.plot = 'all'){ 
  print('reordering...')
  if(proteins.to.plot == 'known'){
    results.dt <- merge(results.dt, metadata.protein[, c('group', 'subgroup', 'Protein.IDs'), with = F], by = 'Protein.IDs', all.x = T)
    results.dt <- results.dt[results.dt$group != 'no group', ]
    protein.ord <- names(sort(table(results.dt$Protein.IDs)))
    pr.order.use.subgroups <- c()
    for (order in kGroupOrder$`Order of appearance`){
      subcomplex <- kGroupOrder[kGroupOrder$`Order of appearance` == order, ]$ComplexName
      proteins.in.subgroup <- unique(results.dt[results.dt$subgroup == subcomplex, ]$Protein.IDs)
      pr.order.use.subgroups <- c(pr.order.use.subgroups, protein.ord[protein.ord %in% proteins.in.subgroup])
    }
    Baits <- rev(c(grep(kTargetProteins[1], protein.ord, value = T),
                   grep(kTargetProteins[2], protein.ord, value = T),
                   grep(kTargetProteins[3], protein.ord, value = T)))
    pr.order.use.subgroups <- unique(c(Baits, pr.order.use.subgroups))
    new.order <- order(factor(results.dt$Protein.IDs, levels = pr.order.use.subgroups, ordered = T))
    results.dt <- results.dt[new.order, ]
    results.dt$Gene <- factor(results.dt$uniprot.gene.best, levels = unique(as.character(results.dt$uniprot.gene.best)))
    group.order <- unique(kGroupOrder[order(kGroupOrder[, c('Order of appearance'), with = F]), ]$Function)
    results.dt$group <- factor(results.dt$group, levels = group.order)
    
  } else if(proteins.to.plot == 'all'){
    target.names <- unique(metadata.protein[metadata.protein$Protein.IDs %in% grep(paste(kTargetProteins, collapse = '|'), 
                                                                                   metadata.protein$Protein.IDs, value = T) , ]$uniprot.gene.best)
    protein.ord <- names(sort(table(results.dt$uniprot.gene.best)))
    protein.ord <- protein.ord[protein.ord %nin% target.names]
    protein.ord <- c(protein.ord, target.names)
    results.dt$Gene <- factor(results.dt$uniprot.gene.best, levels = protein.ord)
  }
  return(results.dt)
}

PlotAmountOfSignificantHeatmap <- function(N.seeds, max.quant.obj, order.experiments.by = 'n.significant', proteins.to.plot = 'all'){
  files <- list.files('./out/QC/pvals_logfc_seeds/', full.names = T)
  metadata.experiment <- max.quant.obj$GetMetadataExperiments()
  metadata.experiment$group <- paste(metadata.experiment$location, metadata.experiment$target, metadata.experiment$buffer, sep = '_')
  metadata.protein <- max.quant.obj$GetMetadataProteins()
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Gene.names
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Protein.IDs
  metadata.protein$uniprot.gene.best <- str_split_fixed(metadata.protein$uniprot.gene.best, ';', 2)[, 1]
  results.list <- list()
  for(file in files){
    seed <- str_split_fixed(file, '_', 5)[, 4]
    dt <- readRDS(file)
    list.of.significant.per.column <- lapply(unique(metadata.experiment$group), function(group) {
      (dt[, paste0("log2.fold.change.", group), with = F] > 1) &
        (dt[, paste0("pvalue.adj.", group), with = F] < 0.05)
    })
    sign.m <- as.data.table(do.call("cbind", list.of.significant.per.column))
    colnames(sign.m) <- unique(metadata.experiment$group)
    sign.m$Protein.IDs <- dt$Protein.IDs
    sign.m <- replace(sign.m, is.na(sign.m), F)
    sign.m <- melt(sign.m, id.vars = 'Protein.IDs', variable.name = 'experiment', value.name = seed)
    results.list[[seed]] <- sign.m
  }
  results.dt <- merge_recurse(results.list)
  results.dt <- results.dt[rowSums(results.dt[, as.character(c(1:N.seeds)), with = F] == TRUE) > 1, ]
  results.dt$n.significant.seeds <- rowSums(results.dt[, as.character(c(1:N.seeds)), with = F] == TRUE)
  results.dt <- results.dt[, c('n.significant.seeds', 'experiment', 'Protein.IDs'), with = F]
  results.dt <- merge(results.dt, metadata.protein[, c('uniprot.gene.best', 'Protein.IDs'), with = F], by = 'Protein.IDs', all.x = T)
  results.dt <- ReorderProteins(results.dt, metadata.protein, proteins.to.plot)
  if(order.experiments.by == 'n.significant'){
    experiment.ord <- names(sort(table(results.dt$experiment)))  
  } else if(order.experiments.by == 'buffer'){
    experiment.ord <- as.character(sapply(kTargets, FUN = function(x) paste0(x, '_', kBuffers)))
    metadata.experiment$experiment.type <- paste(metadata.experiment$location, '_', metadata.experiment$target, sep = '')
    metadata.experiment <- metadata.experiment[, c('experiment.type', 'group'), with = F]
    metadata.experiment <- metadata.experiment[!duplicated(metadata.experiment), ]
    results.dt <- merge(results.dt, metadata.experiment, by.x = 'experiment', by.y = 'group', all.x = T)
    results.dt$experiment.type <- factor(results.dt$experiment.type, kTargets)
  }
  
  results.dt$experiment <- factor(results.dt$experiment, levels = experiment.ord)
  g <- ggplot(results.dt, aes(x = experiment, y = Gene, fill = n.significant.seeds)) +
    geom_tile() +
    scale_fill_viridis(direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major.x = element_blank())
  if(proteins.to.plot == 'known' & order.experiments.by == 'buffer'){
    g <- g +
      facet_grid(rows = c(vars(group)), cols = c(vars(experiment.type)), space = "free", scales = "free")
  } else if(proteins.to.plot == 'known' & !order.experiments.by == 'buffer'){
    g <- g +
      facet_grid(rows = c(vars(group)), space = "free", scales = "free")
  } else if(!proteins.to.plot == 'known' & order.experiments.by == 'buffer'){
    g <- g +
      facet_grid(cols = c(vars(experiment.type)), space = "free", scales = "free")
  }
  pdf(paste0('./out/QC/', MakeDate(), '_amount_of_significant_proteins_all_heatmap.pdf'), width = 20, height = 65)
  print(g)
  dev.off()
}



