source('lib/constants.R')
source('lib/helper_functions.R')

library(ggplot2)
library(ggthemes)
library(reshape2)
library(stringr)
library(viridis)
require(gridExtra)
require(gplots)
library(circlize)
library(gtable)
library(grid)
library(magrittr)

label_wrap_gen <- function(width = 100) {
  function(variable, value) {
    lapply(strwrap(as.character(value), width = width, simplify = FALSE), 
           paste, collapse = "\n")
  }
}

MakeIntensities <- function(max.quant.obj, control.mode) {
  if (control.mode == 'no'){
    intensities <- max.quant.obj$GetIntensityTable()    
  } else if(control.mode == 'logFC'){
    intensities <- max.quant.obj$GetPvaluesLogFC()
    intensities <- as.data.table(intensities)
    intensities <- intensities[, c(grep('log2.fold.change', colnames(intensities), value = T), 
                                   'Protein.IDs'), with = F]
    colnames(intensities) <- gsub('log2.fold.change.', '', colnames(intensities))
    colnames(intensities) <- gsub('.', ' ', colnames(intensities), fixed = T)
  } else if(control.mode == 'subtract'){
    max.quant.obj$
      PowTwo()
    intensities <- max.quant.obj$GetIntensityTable()
    intensities[intensities == 1] <- 0
    metadata.exp <- max.quant.obj$GetMetadataExperiments()
    for (group. in unique(metadata.exp$group)){
      control.exps <- metadata.exp[metadata.exp$group == group. & 
                                     metadata.exp$is.case == F, ]$experiment.name
      case.exps <- metadata.exp[metadata.exp$group == group. & metadata.exp$is.case == T, 
                                ]$experiment.name
      mean.control.exps <- rowMeans(intensities[, control.exps, with = F])
      data.table::set(intensities, j = case.exps, value = intensities[, case.exps, with = F] - 
                        mean.control.exps)
    }
    max.quant.obj$SetIntensityTable(intensities, change.string = 'subtracted background')
    max.quant.obj$LogTricky()
    intensities <- max.quant.obj$GetIntensityTable() 
  } else {
    print('I do not have such control.mode!')
  }
  intensities[intensities <= 0] <- NA
  return(intensities)
}

MakeMean <- function(intensities, metadata.protein, metadata.experiment, control.mode){
  metadata.protein[is.na(metadata.protein$uniprot.gene.best), ]$uniprot.gene.best <- 
    str_split_fixed(metadata.protein[is.na(metadata.protein$uniprot.gene.best), ]$Gene.names, ';', 2)[, 1]
  ids.in.groups <- metadata.protein[metadata.protein$group != 'no group',]$Protein.IDs
  intensities <- intensities[intensities$`Protein IDs` %in% ids.in.groups, ]
  
  intensities[is.na(intensities)] <- 0
  intensities <- merge(intensities, metadata.protein[, c('Protein.IDs', 'uniprot.gene.best', 
                                                         'group', 'subgroup')], 
                       by.x = 'Protein IDs', by.y = 'Protein.IDs')
  
  
  intensities.melt <- melt(intensities, id = c('Protein IDs', 'group', 'subgroup', 'uniprot.gene.best'))
  colnames(metadata.experiment)[colnames(metadata.experiment) == 'group'] <- 'exp' 
  
  if (control.mode == 'logFC') { 
    metadata.experiment <- metadata.experiment[, c("location", "target", 'buffer', 'exp'), with = F]
    metadata.experiment <- metadata.experiment[! duplicated(metadata.experiment), ]
    intensities.melt <- merge(intensities.melt, metadata.experiment, by.x = 'variable', by.y = 'exp', 
                              all.x = T)
    aggr <- intensities.melt
    colnames(aggr) <- c('experiment', 'Protein.IDs', 'group', 'subgroup', 'Gene.name', 'log.Intensity',
                        'location', 'target', 'buffer')
  } else {
    intensities.melt <- merge(intensities.melt, metadata.experiment, by.x = 'variable', 
                              by.y = 'experiment.name', all.x = T)
    intensities.melt <- intensities.melt[intensities.melt$is.case == T, ]
    intensities.melt <- intensities.melt[, names(intensities.melt) %nin% c("variable", "replicate", 
                                                                           "is.case"), with = F]
    aggr <- aggregate(intensities.melt, by = list(intensities.melt$`Protein IDs`, 
                                                  intensities.melt$group, 
                                                  intensities.melt$subgroup, 
                                                  intensities.melt$uniprot.gene.best, 
                                                  intensities.melt$exp,
                                                  intensities.melt$location,
                                                  intensities.melt$buffer,
                                                  intensities.melt$target,
                                                  intensities.melt$pretty_name), FUN = mean, na.rm = T)
    aggr <- aggr[, ! names(aggr) %in% c('Protein IDs', 'group', 'subgroup', 'uniprot.gene.best',
                                        'exp', 'location', 'buffer', 'target', 'pretty_name')]
    colnames(aggr) <- c('Protein.IDs', 'group', 'subgroup', 'Gene.name', 'experiment', 'location', 
                        'buffer', 'target', 'pretty_name', 'log.Intensity')
    aggr$log.Intensity <- aggr$log.Intensity / log2(10)
  }
  if (sum(aggr$log.Intensity == 0) > 0) {
    aggr[aggr$log.Intensity == 0, ]$log.Intensity <- NA 
  }
  return(aggr)
}

ReorderProteinsInHeatmap <- function(aggr, max.quant.obj.tmp, control.mode, metrics, recalculate.dist, 
                            method_ = 'complete'){ 
    print('reordering...')
    if(control.mode == 'subtract' & recalculate.dist == T){
      max.quant.obj.tmp$
        CountDistanceMatrix(metrics = 'eLife', all.targets = T)
    }
    dm <- max.quant.obj.tmp$GetDistanceMatrix()
    dm <- subset(dm, unique(aggr$Protein.IDs))
    hc <- hclust(dm, method = method_)
    #plot(hc)
    pr.order <- dimnames(dm)[hc$order]
    pr.order.use.subgroups <- c()
    
    for (order in kGroupOrder$`Order of appearance`){
      subcomplex <- subset(kGroupOrder, `Order of appearance` == order)$ComplexName
      proteins.in.subgroup <- unique(aggr[aggr$subgroup == subcomplex, ]$Protein.IDs)
      pr.order.use.subgroups <- c(pr.order.use.subgroups, pr.order[pr.order %in% proteins.in.subgroup])
    }
    Baits <- rev(c(grep(kTargetProteins[1], pr.order, value = T),
                   grep(kTargetProteins[2], pr.order, value = T),
                   grep(kTargetProteins[3], pr.order, value = T)))
    pr.order.use.subgroups <- c(Baits, pr.order.use.subgroups)
    new.order <- order(factor(aggr$Protein.IDs, levels = pr.order.use.subgroups, ordered = T))
    aggr <- aggr[new.order, ]
    aggr$Gene.name <- factor(aggr$Gene.name, levels = unique(as.character(aggr$Gene.name)))
    return(aggr)
}

MakeTableOfDistances <- function(aggr){
  aggr.dt <- as.data.frame(aggr)
  aggr.dt[is.na(aggr.dt)] <- 0
  aggr.dt <- as.data.table(aggr.dt)
  aggr.dt$IP <- str_remove(aggr.dt$IP, "_.*")
  groups.list <- unique(aggr$group)
  col.fun = colorRamp2(c(0, 1), c("white", "red"))
  date <- MakeDate()
  pdf(paste0('out/', date, '_distances_on_heatmap.pdf'), width = 8, 15)
  heatmap.list <- lapply(groups.list, function(protein.group) {
    dcasted <- 
      dcast(aggr.dt[aggr.dt$group == protein.group, ], formula = Protein.IDs + buffer ~ IP,
            value.var = "log.Intensity")
    cor.mat <- (1 - cor(dcasted[, as.character(unique(aggr.dt$IP))])) / 2
    melted.cormat <- melt(cor.mat)
    g <- ggplot(data = melted.cormat, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() + 
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name = '(1 - cor)/2') +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) + 
      ggtitle(as.character(protein.group))
    g
  })
  multiplot(plotlist = heatmap.list, layout = matrix(1:10, ncol = 2))
  dev.off()
}

MarkSignificant <- function(aggr){
  significant.dt <- fread('out/significant_proteins.tsv')
  significant.dt <- significant.dt[significant.dt$percentage.significant >= 60, ]
  aggr <- merge(aggr, significant.dt, by = c('Protein.IDs', 'experiment'), all.x = T)
  aggr$is.significant <- T
  aggr[is.na(aggr$percentage.significant), ]$is.significant <- F
  aggr
}


# possible buffers.modes are:
# as is - draw all possible buffers for experiments. Use predefined order 
# add blank - add blank columns if buffer is not present in some experiment
# only shared - deaw only buffers shared by all experiments
# reorder - draw buffers shared by all experiments at the first place

# possible colors are:
# colorful
# BlackandWhite

# possible control.modes are (how to control control or background Intensity):
# no - don't control for background
# subtract - subtract control Intensity from case Intensity
# logFC - log(Intensity_case - Intensity_control)

PlotHeatMap <-
  function(max.quant.obj,
           color = 'colorful',
           buffers.mode = 'as is',
           reorder.by.dist = F,
           make.table.of.distances = F, 
           control.mode = 'no', metrics = 'euclidean', 
           recalculate.dist = F, method_ = "complete") {

  max.quant.obj.tmp <- max.quant.obj$clone()
  metadata <- max.quant.obj.tmp$GetMetadataExperiments()
  
  intensities.no <- MakeIntensities(max.quant.obj.tmp, 'no')
  intensities.sub <- MakeIntensities(max.quant.obj.tmp, 'subtract')
  intensities.log <- MakeIntensities(max.quant.obj.tmp, 'logFC')

  aggr.no <- MakeMean(intensities.no, max.quant.obj.tmp$GetMetadataProteins(), 
                      max.quant.obj.tmp$GetMetadataExperiments(), 'no')
  aggr.sub <- MakeMean(intensities.sub, max.quant.obj.tmp$GetMetadataProteins(), 
                       max.quant.obj.tmp$GetMetadataExperiments(), 'subtract')
  aggr.log <- MakeMean(intensities.log, max.quant.obj.tmp$GetMetadataProteins(), 
                       max.quant.obj.tmp$GetMetadataExperiments(), 'logFC')
  
  min.no <- min(aggr.no$log.Intensity, na.rm = T)
  max.no <- max(aggr.no$log.Intensity, na.rm = T)

  min.sub <- min(aggr.sub$log.Intensity, na.rm = T)
  max.sub <- max(aggr.sub$log.Intensity, na.rm = T)
  
  min.log <- min(aggr.log$log.Intensity, na.rm = T)
  max.log <- max(aggr.log$log.Intensity, na.rm = T) 
  if (control.mode == 'no'){
    aggr <- aggr.no
    min.int <- min(min.no, min.sub)
    max.int <- max(max.no, max.sub)
  } else if (control.mode == 'subtract'){
    aggr <- aggr.sub
    min.int <- min(min.no, min.sub)
    max.int <- max(max.no, max.sub)
  } else if (control.mode == 'logFC'){
    aggr <- aggr.log
    min.int <- min.log
    max.int <- max.log
  } else {
    print('I do not know such control.mode!')
  }
  
  
  aggr$Gene.name <- factor(aggr$Gene.name)
  aggr[grep('Q53F19|P52298|Q09161', aggr$Protein.IDs),]$group <- 'Baits'
  aggr[grep('Q53F19|P52298|Q09161', aggr$Protein.IDs),]$subgroup <- 'Baits'
  aggr <- aggr[aggr$group != 'no group', ]

  group.levels <- c('Baits', unique(kGroupOrder$Function))
  aggr$group <- factor(aggr$group, levels = group.levels)
  
  subgroup.levels <- c('Baits', unique(kGroupOrder$ComplexName))
  aggr$subgroup <- factor(aggr$subgroup, levels = subgroup.levels)
  
  aggr$IP <- metadata[match(aggr$experiment, metadata$group), ]$pretty_name
  buffer.levels <- kBuffers
  aggr$buffer <- factor(aggr$buffer, levels = buffer.levels)
  aggr$IP <- factor(aggr$IP)
  

  if (buffers.mode == 'only shared') {
    ips <- unique(aggr$IP)
    shared.buffers <- unique(aggr[aggr$IP == ips[1], ]$buffer)
    for (ip in ips[2:length(ips)]){
      shared.buffers <- intersect(shared.buffers, unique(aggr[aggr$IP == ip, ]$buffer))
    }
    aggr <- aggr[aggr$buffer %in% shared.buffers, ]
  } else if (buffers.mode == 'reorder') {
    ips <- unique(aggr$IP)
    shared.buffers <- unique(aggr[aggr$IP == ips[1], ]$buffer)
    for (ip in ips[2:length(ips)]){
      shared.buffers <- intersect(shared.buffers, unique(aggr[aggr$IP == ip, ]$buffer))
    }
    buffer.levels <- c(intersect(kBuffers, shared.buffers), setdiff(kBuffers, shared.buffers))
    aggr$buffer <- factor(aggr$buffer, levels = buffer.levels)
  }
  
  if (make.table.of.distances) {
    MakeTableOfDistances(aggr)
  }
  
  if (reorder.by.dist) {
    aggr <- ReorderProteinsInHeatmap(aggr, max.quant.obj.tmp, control.mode, metrics, recalculate.dist, method_)
  }
 
  aggr <- MarkSignificant(aggr)
  pdf(paste0('out/', MakeDate(), '_', color, '_',  buffers.mode, '_', control.mode, 'method=', method_,
             '_heatmap.pdf'), width = 8, height = 30)
  aggr <- data.frame(aggr)
  p <- ggplot(aggr, aes(x = buffer, y = Gene.name, fill = log.Intensity)) +
    ggtitle(paste('Use control: ', control.mode, ', Buffers: ', buffers.mode, sep = '')) + 
    geom_tile(aes(fill = log.Intensity)) +
    geom_point(data = aggr[aggr$is.significant, ], aes(x = buffer, y = Gene.name, colour = 'black'), size = 1, shape = 8) +
    scale_colour_manual(name = '', values = 'black', labels = 'significant')
    
    # geom_raster(aes(fill = log.Intensity)) +
    # geom_tile(aes(colour = is.significant), fill = '#00000000', size = 1) +
    # scale_color_manual(name = "", labels = c('not significant', 'significant'), values = c('#00000000', 'black')) 
   
 
  
  if (color == 'colorful') {
    p <- p + 
      #scale_fill_viridis(option = "A", direction = 1) +
      #scale_fill_gradient(low = "#56B4E9", high = "#E69F00", na.value = "grey50") + 
      scale_fill_gradient(low = "#56B4E9", high = "#D55E00", na.value = "grey50", 
                          limits = c(min.int, max.int)) + 
      theme(panel.background = element_rect(fill = "grey50"))  +
      theme(panel.grid.major = element_blank())
  } else if (color == 'BlackandWhite') {
    p <- p + 
      scale_fill_gradient(low = "light grey", high = "black", na.value = "white", 
                          limits = c(min.int, max.int)) +
      theme_classic()
  } else {
    print('I do not know such color!')
  }
  
  if (buffers.mode == 'add blank') {
    p <- p + 
      facet_grid(rows = vars(group), cols = vars(IP), 
                 space = "free_y", scales = "free_y")
  } else if (buffers.mode == 'as is' | buffers.mode == 'only shared' | buffers.mode == 'reorder') {
    p <- p +
      facet_grid(cols = vars(IP) , rows = c(vars(group)), space = "free", scales = "free")
  } else {
    print('I do not have such buffer.mode!')
  }
  
  if (control.mode == 'logFC'){
    p <- p + 
      guides(fill = guide_colorbar(title = "logFC", reverse = F))
  }
  print(p)
  dev.off()    
  return(0)
  }
