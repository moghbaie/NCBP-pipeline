library(data.table)
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
source('lib/constants.R')
dir.create(file.path("out/", "QC"), showWarnings = FALSE)


PlotIntervalOfIntensities <- function(max.quant.obj, proteins.to.plot){
  metadata <- max.quant.obj$GetMetadataExperiments()
  protein.groups.f <- max.quant.obj$GetIntensityTable()
  protein.groups.f[protein.groups.f == 0] <- NA
  min.intensities <- apply(protein.groups.f, 2, function(column) min(column, na.rm = T))
  min.intensities <- min.intensities[names(min.intensities) != 'Protein IDs']
  min.intensities['Gene.names'] <- 'min intensity'
  max.intensities <- apply(protein.groups.f, 2, function(column) max(column, na.rm = T))
  max.intensities <- max.intensities[names(max.intensities) != 'Protein IDs']
  max.intensities['Gene.names'] <- 'max intensity'
  protein.groups.f <- merge(protein.groups.f,
                            max.quant.obj$GetMetadataProteins()[, c('Protein.IDs', 'Gene.names')],
                            by.x = 'Protein IDs', by.y = 'Protein.IDs')
  protein.groups.f$Gene.names <- str_split_fixed(protein.groups.f$Gene.names, ';', Inf)[, 1]
  protein.groups.f[grep('Q9U6Y5', protein.groups.f$`Protein IDs`), ]$Gene.names <- 'GFP'
  protein.groups.f[protein.groups.f$`Protein IDs` == 'P52298-2', ]$Gene.names <- 'NCBP2.2'
  
  rows.to.plot <- grep(paste(proteins.to.plot, collapse = '|'), protein.groups.f$`Protein IDs`)
  intensities.to.plot <- protein.groups.f[rows.to.plot, ]
  intensities.to.plot <- intensities.to.plot[, -c('Protein IDs')]
  
  intensities.to.plot <- rbind(intensities.to.plot, t(min.intensities))
  intensities.to.plot <- rbind(intensities.to.plot, t(max.intensities))
  intensities.to.plot[is.na(intensities.to.plot)] <- 0
  intensities.to.plot.tmp <- t(intensities.to.plot)
  colnames(intensities.to.plot.tmp) <- intensities.to.plot$Gene.names
  intensities.to.plot <- intensities.to.plot.tmp
  intensities.to.plot.tmp <- NaN
  intensities.to.plot <- intensities.to.plot[-nrow(intensities.to.plot), ]
  intensities.to.plot <- as.data.frame(intensities.to.plot)
  for(col in colnames(intensities.to.plot)){
    intensities.to.plot[[col]] <- as.numeric(as.character(intensities.to.plot[[col]]))
  }
  experiments <- rownames(intensities.to.plot)
  experiments <- c(sort(experiments[-grep('_C_', experiments)]), sort(experiments[grep('_C_', experiments)]))
  intensities.to.plot$experiment <- factor(rownames(intensities.to.plot), levels = experiments)
  
  intensities.to.plot.tmp <- melt(intensities.to.plot, id = c('experiment'))
  colnames(intensities.to.plot.tmp) <- c('experiment', 'protein', 'intensity')
  intensities.to.plot.tmp$intensity <- as.numeric(intensities.to.plot.tmp$intensity)
  
  intensities.to.plot.tmp <- merge(intensities.to.plot.tmp, metadata[, c('experiment.name', 'pretty_name'), with = F], 
                                   by.x = 'experiment', by.y = 'experiment.name')
  intensities.to.plot.tmp$experiment <- gsub('Columbia_NCBP3_', '', intensities.to.plot.tmp$experiment)
  
  intensities.to.plot.tmp$experiment <- gsub('Rockefeller_NCBP[1,2]_', '', intensities.to.plot.tmp$experiment)
  intensities.to.plot.tmp$experiment <- paste(intensities.to.plot.tmp$pretty_name, intensities.to.plot.tmp$experiment, sep = '_')
  
  intensities.to.plot.tmp$protein <- gsub('C17orf85', 'NCBP3', intensities.to.plot.tmp$protein)
  
  
  colours <- unname(kColorList[c(8, 2, 10, 7, 6, 1, 1)])
  names(colours) <- c("NCBP1", "NCBP2", "NCBP2.2", "NCBP3", "GFP", "min intensity", 'max intensity')
  
  dir.create("out/QC/", recursive = T, showWarnings = F)
  
  pdf(paste0("out/QC/",  MakeDate(), '_choose_intensity_to_normalize.pdf'), width = 13, height = 7)
  print(ggplot(intensities.to.plot.tmp, aes(x = experiment, y = intensity, col = protein)) + 
          geom_point() + 
          geom_line(aes(group = protein)) +
          scale_color_manual(values = colours) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          ylab(label = 'log(intensity)') + 
          theme(legend.text = element_text(size = 20)))
  dev.off()
}




