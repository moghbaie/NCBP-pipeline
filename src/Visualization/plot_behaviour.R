require(data.table)
require(stringr)
require(dplyr)

source("src/Clusterization/make_dendrogram.R")
source("lib/helper_functions.R")

PlotBehaviourInner <- function(proteins.list, color.list, df.target, metadata, targets.info, is.first,
                               title = '') {
  for (p in proteins.list) {
    y <- as.numeric(df.target[p, ])
    max.intensity <- max(df.target[proteins.list, ], na.rm = T)
    labels <- str_remove(colnames(df.target), ".*_")
    if (is.first) {
      plot(x = c(1:length(y)), y, ylim = c(0, max.intensity),
           xaxt = "n", type="l", xlab = "",
           ylab = "intensity", bty="n", lwd = 2,
           col = color.list[metadata[Protein.IDs == p, ]$group],
           main = title)
      axis(1, at=c(1:length(y)), labels, las = 2, cex.axis = 0.6)
      for (i in c(1:targets.info$count)) {
        axis(1,at=c(targets.info$lengths[i] + 1, targets.info$lengths[i+1]),
             col="red", line=2.2, tick=T, labels=rep("",2), lwd=2, lwd.ticks=0)
      }
      axis(1, at = targets.info$coords, tick = F, labels = targets.info$names, line=1, cex.axis = 0.5)
      is.first = F
    } else {
      lines(x = c(1:length(y)), y, xaxt = "n", lwd = 2,
            col = color.list[metadata[Protein.IDs == p, ]$group])
    }
  }
  
  return(is.first)
}

PlotBehaviour <- function(proteins.list, color.list, df.target, metadata, title = '', targets.info) {
  # cluster list contains vectors of protein names
  
  is.first = T
  metadata.c <- metadata[metadata$Protein.IDs %in% proteins.list, ]
  grey.proteins <- metadata.c[group == "no group", ]$Protein.IDs
  known.proteins <- metadata.c[group != "no group", ]$Protein.IDs
    
  is.first <- 
    PlotBehaviourInner(grey.proteins, color.list, df.target, metadata, targets.info, is.first, title)
  PlotBehaviourInner(known.proteins, color.list, df.target, metadata, targets.info, is.first, title)
}

CountTargetsLengthsAndCoords <- function(max.quant.obj) {
  targets <- unique(paste(max.quant.obj$GetMetadataExperiments()$location, 
                          max.quant.obj$GetMetadataExperiments()$target, sep = "_"))
  
  target.cases <- list()
  for (current.target in targets) {
    target.cases[[current.target]] <- 
      max.quant.obj$GetMetadataExperiments()[intersect(grep(current.target, experiment.name),
                                                    which(is.case)), ]$experiment.name
  }
  
  lengths.targets <- as.numeric(cumsum(lapply(target.cases, length)))
  lengths.targets <- c(0, lengths.targets)
  target.coords <- c()
  for (i in c(1:length(targets))) {
    target.coords <- c(target.coords, (lengths.targets[i] + lengths.targets[i + 1])/2 + 1)
  }
  
  return(list(lengths = lengths.targets, coords = target.coords, 
              count = length(targets), names = targets))
}

PlotDendrosWithBehaviour <- function(max.quant.obj, n.of.clusters = 9) {

  clusters.target <- max.quant.obj$GetClusters()
  df.intensities <- MakeDataFrameWithCases(max.quant.obj$GetIntensityTable(), 
                                           id.col = kIdColumnName, 
                                           metadata = max.quant.obj$GetMetadataExperiments(),
                                           keep.id.column = F)
  
  df.intensities[is.na(df.intensities)] <- 0
  count.clusters <- length(unique(clusters.target))
  
  targets.info <- CountTargetsLengthsAndCoords(max.quant.obj)
  metadata.experiments <- max.quant.obj$GetMetadataExperiments()
  metadata.experiments$exp <- paste(metadata.experiments$location, metadata.experiments$target, sep = '_')
  targets.info$names <- metadata.experiments[match(targets.info$names, metadata.experiments$exp)]$pretty_name
  
  pdf(paste0("out/", MakeDate(), "_behaviour_plots.pdf"), width = 10, height = 20)
  layout(matrix(c(1, 1, 2, 3), 2, byrow = T),
         widths = c(.85, .15), heights = c(.85, .15))
  
  list.proteins.in.clusters <- list()
  
  for (cluster in c(1:count.clusters)) {
    proteins.in.cluster <- names(clusters.target[clusters.target == cluster])
    list.proteins.in.clusters[[cluster]] <- proteins.in.cluster
    if (length(proteins.in.cluster) < 2) {
      plot.new()
      par(mar = c(4,2,4,0))
      
      PlotBehaviour(proteins.list = proteins.in.cluster, color.list = kColorList,
                    df.target = df.intensities, metadata = max.quant.obj$GetMetadataProteins(),
                    title = cluster, targets.info = targets.info)
      
      par(mar = c(2,0,3,0))
      plot.new()
      legend(x = 0, y = 1, col = kColorList, cex = 0.8,
             legend = names(kColorList), bty = 'n', pch = 20)
      next
    }
    
    dist.mtx.cluster <- 
      as.dist(as.matrix(max.quant.obj$GetDistanceMatrix())[proteins.in.cluster, 
                                                           proteins.in.cluster])
    
    dendro.one.cluster <- 
      MakeDendroWithClusters(dist.mtx.cluster, max.quant.obj$GetMetadataProteins(), 
                             n.of.clusters = 1)
    
    par(mar = c(4,1,1,4))
    
    plot(dendro.one.cluster$dend, horiz=T)
    
    par(mar = c(4,2,4,0))
    
    PlotBehaviour(proteins.list = proteins.in.cluster, color.list = kColorList,
                  df.target = df.intensities, metadata = max.quant.obj$GetMetadataProteins(),
                  title = cluster, targets.info = targets.info)
    
    par(mar = c(2,0,3,0))
    plot.new()
    legend(x = 0, y = 1, col = kColorList, cex = 0.8,
           legend = names(kColorList), bty = 'n', pch = 20)
  }
  dev.off()
}
