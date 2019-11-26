library(gridExtra)
library(grid)
library(ggplot2)

PlotWave <- function(max.quant.obj, color.list = kColorList, n.anovas.passed = c(1, 2)) {
  metadata <- max.quant.obj$GetMetadataProteins()
  metadata.exp <- max.quant.obj$GetMetadataExperiments()
  metadata.exp <- metadata.exp[, c("group", "buffer", "location", "target"), with = F]
  metadata.exp <- metadata.exp[!duplicated(metadata.exp), ]
  metadata.experiment <- max.quant.obj$GetMetadataExperiments()
  
  max.quant.obj <- max.quant.obj$MakeObjectOnlySignificant()
  
  pvals.logfc.1 <- max.quant.obj$GetPvaluesLogFC()
  pvals.logfc.1 <- pvals.logfc.1[pvals.logfc.1$n.of.anovas.passed >= n.anovas.passed[1], ]
  logfc.cols.1 <- colnames(pvals.logfc.1)[grep("log", colnames(pvals.logfc.1))]
  pval.cols.1 <- colnames(pvals.logfc.1)[grep("pval", colnames(pvals.logfc.1))]
  
  if (length(n.anovas.passed) > 1) {
    pvals.logfc.2 <- max.quant.obj$GetPvaluesLogFC()
    pvals.logfc.2 <- pvals.logfc.2[pvals.logfc.2$n.of.anovas.passed >= n.anovas.passed[2], ]
    logfc.cols.2 <- colnames(pvals.logfc.2)[grep("log", colnames(pvals.logfc.2))]
    pval.cols.2 <- colnames(pvals.logfc.2)[grep("pval", colnames(pvals.logfc.2))]
  }
  
  
  targets <- sort(unique(max.quant.obj$GetMetadataExperiments()$target))
  
  protein.groups <- unique(metadata[group != "no group", group])
  
  n.of.groups <- length(protein.groups)
  experiment.names <- unique(max.quant.obj$GetMetadataExperiments()$group)
  
  pdf(paste0('out/', MakeDate(), "_signif_proteins_in_complexes.pdf"), width = 7, height = 20)
  plot.data <- list()
  
  for (group.name in protein.groups) {
    dt <- rbindlist(lapply(experiment.names, function(experiment) {
      
      group.protein.ids.1 = unlist(metadata[(group == group.name), "Protein.IDs", with = F])
      pvals.logfc.g.1 <- as.data.table(pvals.logfc.1[pvals.logfc.1$Protein.IDs %in% group.protein.ids.1, ])
      logfc.col.1 <- logfc.cols.1[grep(paste0(experiment, "$"), logfc.cols.1)]
      pval.col.1 <- pval.cols.1[grep(paste0(experiment, "$"), pval.cols.1)]
      n_of_significant.1 <- sum((pvals.logfc.g.1[, pval.col.1, with = F] < 0.05) & 
                                (pvals.logfc.g.1[, logfc.col.1, with = F] > 1), na.rm = T)
      percent_of_significant.1 <- n_of_significant.1 / nrow(metadata[(group == group.name), ])
      
      if(length(n.anovas.passed) > 1){
        group.protein.ids.2 = unlist(metadata[(group == group.name), "Protein.IDs", with = F])
        pvals.logfc.g.2 <- as.data.table(pvals.logfc.2[pvals.logfc.2$Protein.IDs %in% group.protein.ids.2, ])
        logfc.col.2 <- logfc.cols.2[grep(paste0(experiment, "$"), logfc.cols.2)]
        pval.col.2 <- pval.cols.2[grep(paste0(experiment, "$"), pval.cols.2)]
        n_of_significant.2 <- sum((pvals.logfc.g.2[, pval.col.2, with = F] < 0.05) & 
                                  (pvals.logfc.g.2[, logfc.col.2, with = F] > 1), na.rm = T)
        percent_of_significant.2 <- n_of_significant.2 / nrow(metadata[(group == group.name), ])
        data.table(experiment = experiment, n_of_significant.1 = n_of_significant.1,
                   percent_of_significant.1 = percent_of_significant.1, 
                   n_of_significant.2 = n_of_significant.2,
                   percent_of_significant.2 = percent_of_significant.2)
      } else {
        data.table(experiment = experiment, n_of_significant.1 = n_of_significant.1,
                   percent_of_significant.1 = percent_of_significant.1)
      }
      
    }))
    
    plot.data[[group.name]] <- dt
  }
  
  order.table <- rbindlist(lapply(protein.groups, function(group.name) {
    dt <- plot.data[[group.name]]
    
    normalized.sums <- c()
    for (target in targets) {
      normalized.sums <- c(normalized.sums, sum(dt[grep(target, experiment), percent_of_significant.1])/
                            length(dt[grep(target, experiment), percent_of_significant.1]))
    }
    names(normalized.sums) <- targets
    
    data.table(group = group.name, 
               max.target = names(which.max(normalized.sums)), max.value = max(normalized.sums))
  }))
  # sort complexes
  sorted.protein.groups <- c()
  for (target in targets) {
    dt <- order.table[grep(target, max.target), ]
    setorder(dt, -max.value)
    sorted.protein.groups <- c(sorted.protein.groups, unlist(dt$group))
  }
  layout(matrix(c(1:n.of.groups, rep(n.of.groups + 1, n.of.groups)), ncol = 2, byrow = F),
         widths = c(.83, .17))
  par(mar = c(8, 4, 5, 0)) 
  for (group.name in sorted.protein.groups) { 
    dt <- plot.data[[group.name]]
    # sort buffers
    dt <- merge(dt, metadata.exp, by.x = "experiment", by.y = 'group')
    dt <- dt[order(dt$target, dt$location, match(dt$buffer, kBuffers)), ]
    dt$experiment_pretty <- metadata.experiment[match(dt$experiment, metadata.experiment$group), ]$pretty_name
    dt$experiment_pretty <- paste(dt$experiment_pretty, dt$buffer, sep = '_')
    plot(x = c(1:dim(dt)[[1]]), dt$percent_of_significant.1, 
           ylim = c(0, 1),
           xaxt = "n", type="l", xlab = "",
           ylab = "Proportion of significant proteins", bty="n", lwd = 2,
           col = color.list[group.name],
           main = paste("Proportion of significant proteins from group", group.name))
    if(length(n.anovas.passed) > 1){
      lines(x = c(1:dim(dt)[[1]]), dt$percent_of_significant.2, lty = 2,
            col = color.list[group.name])
    }
  
    axis(1, at = c(1:dim(dt)[[1]]), 
         labels = dt$experiment_pretty, las = 2, cex.axis = 0.6)
    # par(new = T)
    # plot(x = c(1:dim(dt)[[1]]), dt$n_of_significant.1, 
    #      axes = F, xlab = NA, ylab = NA, col = color.list[group.name], pch = 20, cex = 0.01)
    # axis(side = 4)
    # mtext(side = 4, 'Number of significan proteins', cex = 0.7, line = 3)
    
  }
  par(mar = c(1, 0, 5, 0))
  plot.new()
  # legend("topleft", col = color.list, cex = 1,
  #        legend = names(color.list), bty = 'n', pch = 20)
  # if(length(n.anovas.passed) > 1){
  #   legend("topleft", inset = .1, col = "black", cex = 1,
  #          legend = c("once", "twice"), lty = c(1, 2), 
  #          lwd = c(2, 2), title = "passed ANOVA", bty = 'n')
  # }
  metadata.to.print <- metadata[metadata$group != 'no group', ]
  metadata.to.print[is.na(metadata.to.print$uniprot.gene.best), ]$uniprot.gene.best <- 
    str_split_fixed(metadata.to.print[is.na(metadata.to.print$uniprot.gene.best), ]$Gene.names, ';', 2)[, 1]
  metadata.to.print <- metadata.to.print[, c('uniprot.gene.best', 'group'), with = F]
  colnames(metadata.to.print) <- c('gene name', 'group')
  Xs <- rep(0, length(unique(metadata.to.print$group)))
  Ys <- seq(1, 0, length.out = length(unique(metadata.to.print$group)))
  cols <- color.list[match(sorted.protein.groups, names(color.list))]
  labels. <- 
    unlist(lapply(sorted.protein.groups, function(group.) {
      if(group. == 'Spliceosome'){
        x <- metadata.to.print[metadata.to.print$group == group., ]$`gene name`
        y <- c()
        for(i in c(1:length(x))){
          y <- c(y, x[i], ',')
          if(i %% 4 == 0){
            y <- c(y, '\n')
          }
        }
        y <- paste(y, collapse = '')
        substr(y, 1, (nchar(y) - 1))
      } else{
        paste(metadata.to.print[metadata.to.print$group == group., ]$`gene name`, collapse = '\n')   
      }
      
    }))
  
  text(x = Xs, y = Ys,   
       labels = labels., 
       col = cols,                 
       cex = 0.5, adj = 0)          
  dev.off()
  
 }

