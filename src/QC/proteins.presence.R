library(data.table)
library(stringr)
library(bit64)
library(dplyr)
library(ggplot2)
library(reshape2)

dir.create(file.path("out/", "QC"), showWarnings = FALSE)

PlotTargetIntensityHist <- function(max.quant.obj, proteins.to.plot, new.data = F){
  intensities <- max.quant.obj$GetIntensityTable()
  metadata.exp <- max.quant.obj$GetMetadataExperiments()
  metadata.prot <- max.quant.obj$GetMetadataProteins()
  intensities <- merge(intensities, metadata.prot[, c("Protein.IDs", 'Gene.names')], 
                       by.x = 'Protein IDs', by.y = "Protein.IDs")
  intensities$Gene.names <- str_split_fixed(intensities$Gene.names, ';', Inf)[, 1]
  intensities[grep('Q53F19', intensities$`Protein IDs`), ]$Gene.names <- 'NCBP3'
  pdf(paste("out/QC/", MakeDate(), '_targets_int_density.pdf'), width = 7, height = 5)
  par(mfrow = c(2, 2))
  for (protein in proteins.to.plot){
    protein.name <- intensities[grep(protein, intensities$`Protein IDs`), ]$Gene.names[1]
    print(protein.name)
    if (protein.name == 'NCBP3'){
      if(new.data){
        location.s <- c("Columbia")
      } else{
        location.s <- c('Rockefeller', "Columbia")
      }
    } else{ 
      location.s <- 'Rockefeller'
    }
    for (location. in location.s){
      print(location.)
      intensities.tmp <- intensities[, grep('_C_', grep(location., grep(protein.name, 
                                                                       colnames(intensities), value = T), value = T), 
                                            value = T, invert = T), with = F]
      intensities.tmp2 <- intensities[, c(grep('_C_', grep(location., grep(protein.name, 
                                                                       colnames(intensities), value = T), value = T), 
                                            value = T, invert = T), 'Protein IDs'), with = F]
      intensities.tmp[intensities.tmp == 0] <- NA
      intensities.tmp <- intensities.tmp[!rowSums(is.na(intensities.tmp)) == ncol(intensities.tmp), ]
      dens <- apply(intensities.tmp, 2, density, na.rm = T)
      target.intensities <- as.numeric(unname(unlist(intensities.tmp2[grep(paste0(protein, '\\;|', protein, '$'), 
                                                                           intensities.tmp2$`Protein IDs`, value = T), ])))
      target.intensities <- target.intensities[! is.na(target.intensities)]
      means <- colMeans(intensities.tmp, na.rm = T)[target.intensities != 0]
      sds <- unlist(apply(intensities.tmp, 2, function(col) sd(col, na.rm = T)))[target.intensities != 0]
      target.intensities <- target.intensities[target.intensities != 0]
      delta.means <- round(mean(abs(means - target.intensities) / sds), 1)
      plot(NA, xlim = range(sapply(dens, "[", "x")), ylim = range(sapply(dens, "[", "y")), 
           main = unique(metadata.exp[metadata.exp$location == location. & metadata.exp$target == protein.name, ]$pretty_name), 
           xlab = 'Log Intensity', ylab = 'Density', sub = paste('z-score=', delta.means))
      mapply(lines, dens, col = "blue")
      abline(v = target.intensities, col = "red")       
    }
  }
  dev.off()
}
