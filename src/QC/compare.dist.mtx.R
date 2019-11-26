library(ade4)
library(stringr)
library(cba)
dir.create(file.path("out/", "QC"), showWarnings = FALSE)
CompareDistanceMtxs <- function(dist.mtx1, dist.mtx2){
  common.names <- intersect(names(dist.mtx1), names(dist.mtx2))
  dist.mtx1 <- subset(dist.mtx1, which(names(dist.mtx1) %in% common.names))
  dist.mtx2 <- subset(dist.mtx2, which(names(dist.mtx2) %in% common.names))
  r1 <- mantel.randtest(dist.mtx1, dist.mtx2)
  p <- r1$plot
  r1$plot$xlim <- c(min(p$xlim), (max(p$xlim) + 0.1))
  pdf(paste('out/QC/', MakeDate(), '_distance_mtx_similarity.pdf', sep = ''), height = 5, width = 7)
  plot(r1, main = "Mantel's test", xlab = 'simiarity')
  mtext(paste('p-value = ', r1$pvalue, sep = ''), side = 3)
  dev.off()
  return(r1)
}


