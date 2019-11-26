library(stringr)
library(data.table)
library(bit64)
library(dplyr)
library(dendextend)
library(gplots)
library(ggplot2)


dir.create(file.path("out/", "QC"), showWarnings = FALSE)

PlotTICvsContaminants <- function(path.to.mq){
  metadata.exp <- fread('data/pretty_experiments_names.csv')
  msscans <- fread(paste(path.to.mq, "msmsScans.txt", sep = ''))
  msscans <- msscans[, c("Experiment", "Scan number", "Total ion current")]
  TIC <- aggregate(`Total ion current` ~ `Experiment`, msscans, sum)
  TIC$Experiment <- gsub('tk', 'Columbia', TIC$Experiment)
  TIC$Experiment <- gsub('^NCBP', 'Rockefeller_NCBP', TIC$Experiment)
  prot.groups <- fread(paste(path.to.mq, "proteinGroups.txt", sep = ''))
  metadata <- prot.groups[, c(1:12), with = F]
  metadata[, potential.contaminant := prot.groups$`Potential contaminant`]
  intensity_cols <- colnames(prot.groups[, grep("iBAQ *.",
                                                colnames(prot.groups)), with=F])
  intensities <- prot.groups[, intensity_cols, with = F] / 100000
  names(intensities) <- str_remove(names(intensities), "iBAQ")
  colnames(intensities) <- gsub('^ tk', 'Columbia', colnames(intensities))
  colnames(intensities) <- gsub('^ NCBP', 'Rockefeller_NCBP', colnames(intensities))
  metadata[, contaminant.bool := ifelse(potential.contaminant == "+", T, F)]
  sum.intensities <- colSums(intensities, na.rm = T)
  sum.int.contam <- colSums(intensities[metadata$contaminant.bool, ], na.rm = T)
  contam.abundance <- sum.int.contam / sum.intensities
  TIC <- TIC[order(TIC$Experiment), ]
  contam.abundance <- contam.abundance[order(names(contam.abundance))]
  contam.TIC <- data.frame(contam.abundance = unname(contam.abundance), TIC = TIC$`Total ion current`, protein = TIC$Experiment)
  
  date <- MakeDate()
  pdf(paste("out/QC/", date, "_TIC_vs_contam_abund.pdf", sep = ''), height = 6, width = 10)
  contam.TIC$TIC <- log(contam.TIC$TIC)
  plot(x = contam.TIC$contam.abundance, y = contam.TIC$TIC, xlab = "Contaminants abundance", ylab = "log(Total ion current)")
  exps <- unique(paste(str_split_fixed(contam.TIC$protein, '_', Inf)[, 1], 
                           str_split_fixed(contam.TIC$protein, '_', Inf)[, 2], sep = "_"))
  # CTRL
  points(contam.TIC[grep("_C_", contam.TIC$protein), ]$contam.abundance,
         contam.TIC[grep("_C_", contam.TIC$protein), ]$TIC,
         pch = 22, cex = 1.6, col = "black")
  i = 1
  cors <- c()
  for (exp in exps){
    points(contam.TIC[grep(exp, contam.TIC$protein), ]$contam.abundance,
           contam.TIC[grep(exp, contam.TIC$protein), ]$TIC,
           pch = 20, cex = 1.6, col = kColExp[i])
    i = i + 2
    cors <- c(cors, round(cor(contam.TIC[grep(exp, contam.TIC$protein), ]$contam.abundance,
                        contam.TIC[grep(exp, contam.TIC$protein), ]$TIC), digits = 2))
  }
  
  exps <- unlist(lapply(exps, function(x) metadata.exp[metadata.exp$experiment == x, ]$pretty_name))
  legend("topright", legend = paste(exps, "cor =", cors, sep = ' '), 
         col = kColExp[c(1, 3, 5, 7)], pch = 19, bty = "n")
  dev.off()  
  return(0)
}





