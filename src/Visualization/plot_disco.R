library(data.table)
library(stringr)
library(bit64)
library(dplyr)
library(dendextend)
library(ggplot2)

source('lib/constants.R')
source("lib/helper_functions.R")

# plot 2 histograms
plotTwoExpDT <- function(dt, x1, x2, x1.name, x2.name, 
                         highlite = NULL, highlite.col = NULL,
                         highlite.names = NULL,
                         highlite.legend = NULL, col.lne = NULL, 
                         color.list. = kColorList){
  # default colors
  col.pnt <- rgb(0, 0, 0, 1)
  col.hst <- "white"
  col.hst.sig <- "gray95"
  col.hst.brd <- "black"
  if (is.null(col.lne)) {
    col.lne <- "gray50" 
  }
  orf1.col <- "#D55E00"
  orf2.col <- "#0072B2"
  if (is.null(highlite.col)) {
    highlite.col <- rep('blueviolet', length(highlite))
  }
  
  # pg of ORF1 and ORF2 proteins
  orf1.orf2 <- c(orf1 = "57", 
                 orf2 = "59")
  
  # create histogram objects
  hst.bns <- seq(min(unlist(dt[which(dt[, x1, with = F] != 0), x1, with = F]), 
                     unlist(dt[which(dt[, x2, with = F] != 0), x2, with = F])), 
                 max(unlist(dt[which(dt[, x1, with = F] != 0), x1, with = F]), 
                     unlist(dt[which(dt[, x2, with = F] != 0), x2, with = F])), 
                 length.out = 51)
  h1 <-  hist(unlist(dt[which(dt[, x1, with = F] != 0), x1, with = F]), plot = F, breaks = hst.bns)
  h2 <-  hist(unlist(dt[which(dt[, x2, with = F] != 0), x2, with = F]), plot = F, breaks = hst.bns)
  
  # prepare layout for figure
  layMat <- matrix(c(1, 3, 2, 4, 4, 4), ncol = 2, byrow = F)
  layout(layMat, heights=c(2/8, 4/8, 2/8), widths = c(0.9, 0.1))
  
  # plot histograms (as barplots)
  par(mar = c(.4, 4, 3, 1))
  hist(unlist(dt[which(dt[, x1, with = F] != 0), x1, with = F]), 
       breaks = hst.bns, main = "", las = 1,
       ylim = c(0, round(max(c(h1$counts, h2$counts) + 5), - 1)),
       col = col.hst.sig, xaxt = "n")
  par(xpd = T)
  text(x = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       y = 0, pos = 1, cex = .7, offset = .2,
       c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
       col = "gray20")
  par(xpd = F)
  mtext(sprintf("%s", x1.name), side = 3, line = 0.5, cex = 1)
  
  par(mar=c(3, 4, .4, 1))
  hist(unlist(dt[which(dt[, x2, with = F] != 0), x2, with = F]), 
       breaks = hst.bns, main = "", las = 1, xaxt = "n",
       ylim = c(round(max(c(h1$counts, h2$counts) + 5), -1), 0),
       col = col.hst.sig)
  par(xpd = T)
  text(x = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 
       y = 0, pos = 3, cex = .7, offset = .2,
       c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
       col = "gray20")
  par(xpd = F)
  mtext(x2.name, side = 1, line = 0.5, cex = 1)
  
  # plot middle figure with protein links
  par(mar = c(0, 4, 0, 1))
  plot(NA, type = "n", xlim = c(0, max(unlist(dt[which(dt[, x1, with = F] != 0), x1, with = F]), 
                                       unlist(dt[which(dt[, x2, with = F] != 0), x2, with = F]))), 
       ylim = c(0, 1), axes = F, ann = F)
  
  for (i in 1: nrow(dt)) {
    if (col.lne[i] == "#BEBEBE4D") {
      lines(x = unlist(dt[i, c(x1, x2), with = F]), 
            y = c(1, 0), col = col.lne[i], lwd = 1.5) 
    }
  }
  
  for (i in 1:nrow(dt)) {
    if (col.lne[i] != "#BEBEBE4D") {
      lines(x = unlist(dt[i, c(x1, x2), with = F]), 
            y = c(1, 0), col = col.lne[i], lwd = 1.5) 
    }
  }
  
  text.line <- 0
  text.line.y <- seq(0.9, 0.1, by = - 0.2)
  text.mline <- length(text.line.y)
  
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("center", col = unname(color.list.[c(1: 4)]), cex = 1.5,
         legend = c('not signficant', 'significant in \n upper plot', 
                    'significant in \n lower plot', 'significant in \n both'),
         bty = 'n', pch = 20, y.intersp = 2)
}


PlotDiscoPlots <- function(max.quant.obj, color.list. = kColorList, mode = 'buffer'){
  dt <- max.quant.obj$GetIntensityTable()
  dt[is.na(dt)] <- 0
  dt.mean <- as.data.table(CountMeanBetweenReplicasForCases(dt, 
                                               id.col = "Protein IDs", 
                                               metadata = max.quant.obj$GetMetadataExperiments(),
                                               keep.id.column = T))

  pdf(paste0("out/", MakeDate(), "_mirror_plots_", mode, ".pdf"), width = 15, height = 7)  
 
  if(mode == 'buffer'){
    items <- kBuffers
  } else if (mode == 'target'){
    items <- kTargets
  }

  for (item in items) {
    print(item)
    cols.item <- colnames(dt.mean)[grep(item, colnames(dt.mean))]
    if (length(cols.item) > 2) {
      comb.item <- combn(cols.item, 2)
    } else if (length(cols.item) == 2) {
      comb.item <- as.matrix(cols.item)
    } else {
      next
    }
    for (i in c(1: dim(comb.item)[[2]])) {
      pvals.tmp <-
        as.data.table(max.quant.obj$GetPvaluesLogFC())[, grep(paste(comb.item[1, i], '|', comb.item[2, i], '|Protein.IDs', sep = ''), 
                                                              colnames(max.quant.obj$GetPvaluesLogFC()), value = T), with = F]
      dt.mean.tmp <- dt.mean
      sign.in.x1 <-
        pvals.tmp[which((pvals.tmp[, paste('pvalue.adj_', comb.item[1, i], sep = ''), with = F] < 0.05) &
                          (pvals.tmp[, paste('log2.fold.change_', comb.item[1, i], sep = ''), with = F] > 1)), ]$Protein.IDs
      sign.in.x2 <-
        pvals.tmp[which((pvals.tmp[, paste('pvalue.adj_', comb.item[2, i], sep = ''), with = F] < 0.05) &
                          (pvals.tmp[, paste('log2.fold.change_', comb.item[2, i], sep = ''), with = F] > 1)), ]$Protein.IDs
      sign.in.both <- intersect(sign.in.x1, sign.in.x2)
      dt.mean.tmp$type <- color.list.[1]
      dt.mean.tmp[dt.mean.tmp$Protein.IDs %in% sign.in.x1,]$type <- color.list.[2]
      dt.mean.tmp[dt.mean.tmp$Protein.IDs %in% sign.in.x2,]$type <- color.list.[3]
      dt.mean.tmp[dt.mean.tmp$Protein.IDs %in% sign.in.both,]$type <- color.list.[4]
      line.colors <- dt.mean.tmp$type
      metadata <- max.quant.obj$GetMetadataExperiments()
      
      name1 <- paste(unique(metadata[metadata$group == comb.item[1, i], ]$pretty_name), 
                     unique(metadata[metadata$group == comb.item[1, i], ]$buffer), sep = '_') 
      name2 <- paste(unique(metadata[metadata$group == comb.item[2, i], ]$pretty_name), 
                     unique(metadata[metadata$group == comb.item[2, i], ]$buffer), sep = '_') 
      plotTwoExpDT(dt.mean, comb.item[1, i], comb.item[2, i], 
                   name1, name2, col.lne = line.colors, 
                   color.list. = color.list.)
    }
  }  
 
  dev.off()  
}



