library(reshape)
MaxQuant$set("public", "AnovaVote", function(sign.cutoff = 60, N.seeds) {
  files <- list.files('./out/QC/pvals_logfc_seeds/', full.names = T)
  metadata.experiment <- private$metadata.experiments
  metadata.experiment$group <- paste(metadata.experiment$location, metadata.experiment$target, metadata.experiment$buffer, sep = '_')
  metadata.protein <- private$metadata.proteins
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Gene.names
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Protein.IDs
  metadata.protein$uniprot.gene.best <- str_split_fixed(metadata.protein$uniprot.gene.best, ';', 2)[, 1]
  significant.list <- list()
  pvalue.logfc.list <- list()
  for(file in files){
    seed <- str_split_fixed(file, '_', 5)[, 4]
    dt <- readRDS(file)
    pvalue.logfc <- melt(dt, id.vars = c('Protein.IDs', 'GENES'), value.name = seed)
    pvalue.logfc.list[[seed]] <- pvalue.logfc
    
    list.of.significant.per.column <- lapply(unique(metadata.experiment$group), function(group) {
      (dt[, paste0("log2.fold.change.", group), with = F] > 1) &
        (dt[, paste0("pvalue.adj.", group), with = F] < 0.05)
    })
    sign.m <- as.data.table(do.call("cbind", list.of.significant.per.column))
    colnames(sign.m) <- unique(metadata.experiment$group)
    sign.m$Protein.IDs <- dt$Protein.IDs
    sign.m <- replace(sign.m, is.na(sign.m), F)
    sign.m <- melt(sign.m, id.vars = 'Protein.IDs', variable.name = 'experiment', value.name = seed)
    significant.list[[seed]] <- sign.m
  }
  pvalue.logfc.dt <- merge_recurse(pvalue.logfc.list)
  pvalue.logfc.dt$mean <- unlist(apply(pvalue.logfc.dt[, as.character(c(1:N.seeds)), with = F], 1, function(row) mean(row, na.rm = T)))
  pvalue.logfc.dt$sd <- unlist(apply(pvalue.logfc.dt[, as.character(c(1:N.seeds)), with = F], 1, function(row) sd(row, na.rm = T)))
  pvalue.logfc.dt <- pvalue.logfc.dt[!is.na(pvalue.logfc.dt$sd), ]
  pvalue.logfc.dt <- pvalue.logfc.dt[, c('Protein.IDs', 'GENES', 'variable', 'mean', 'sd'), with = F]
  pvalue.logfc.dt$stat <- 'log2.fold.change'
  pvalue.logfc.dt[grepl('pvalue.adj', pvalue.logfc.dt$variable), ]$stat <- 'pvalue.adj'
  pvalue.logfc.dt$experiment <- gsub('log2.fold.change.', '', pvalue.logfc.dt$variable)
  pvalue.logfc.dt$experiment <- gsub('pvalue.adj.', '', pvalue.logfc.dt$experiment)
  pvalue.logfc.dt <- pvalue.logfc.dt[, c("Protein.IDs", "GENES", "mean", "sd", "stat", "experiment" ), with = F]
  write.table(pvalue.logfc.dt, 'out/pvalue_logfc.tsv', sep = '\t', row.names = F, quote = F)
  metadata.experiment$experiment <- unlist(lapply(metadata.experiment$experiment.name, function(x) substr(x, 1, (nchar(x) - 1))))
  pvalue.logfc.dt$pretty_name <- metadata.experiment[match(pvalue.logfc.dt$experiment, metadata.experiment$experiment), ]$pretty_name
  pvalue.logfc.dt$experiment_pretty <- paste(pvalue.logfc.dt$pretty_name, str_split_fixed(pvalue.logfc.dt$experiment, '_', 3)[, 3], sep = '_')
  
  pvalue.logfc.dt.logFC <- pvalue.logfc.dt[pvalue.logfc.dt$stat == 'log2.fold.change', 
                                           c("Protein.IDs", "GENES", "mean", "sd", "experiment_pretty"), with = F]
  pvalue.logfc.dt.pval <- pvalue.logfc.dt[pvalue.logfc.dt$stat == 'pvalue.adj', 
                                           c("Protein.IDs", "GENES", "mean", "sd", "experiment_pretty"), with = F]
  colnames(pvalue.logfc.dt.logFC)[colnames(pvalue.logfc.dt.logFC) == 'mean'] <- 'log2FC_mean'
  colnames(pvalue.logfc.dt.logFC)[colnames(pvalue.logfc.dt.logFC) == 'sd'] <- 'log2FC_sd'
  colnames(pvalue.logfc.dt.pval)[colnames(pvalue.logfc.dt.pval) == 'mean'] <- 'pval.adj_mean'
  colnames(pvalue.logfc.dt.pval)[colnames(pvalue.logfc.dt.pval) == 'sd'] <- 'pval.adj_sd'
  pvalue.logfc.dt.tmp <- merge(pvalue.logfc.dt.logFC, pvalue.logfc.dt.pval)
  write.table(pvalue.logfc.dt.tmp, 'out/pvalue_logfc_all.tsv', sep = '\t', row.names = F, quote = F)
  
  pvalue.logfc.dt <- dcast(pvalue.logfc.dt, Protein.IDs + GENES ~ stat + experiment, value.var = 'mean')
  private$pvals.logfc <- as.data.table(pvalue.logfc.dt)
  
  significant.dt <- merge_recurse(significant.list)
  significant.dt$percentage.significant <- rowSums(significant.dt[, as.character(c(1:N.seeds)), with = F] == TRUE) / N.seeds * 100
  significant.dt <- significant.dt[, c('Protein.IDs', 'experiment', 'percentage.significant'), with = F]  
  significant.dt <- significant.dt[significant.dt$percentage.significant >= sign.cutoff, ]
  write.table(significant.dt, 'out/significant_proteins.tsv', sep = '\t', row.names = F, quote = F)
  
  private$AddChange("ANOVAs voted")
})
