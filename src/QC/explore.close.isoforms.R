library(ggplot2)
dir.create(file.path("out/", "QC"), showWarnings = FALSE)

ExploreCloseIsoforms <- function(path.to.mq, isoform1, isoform2, old.data = F){
  msscans <- fread(paste(path.to.mq, "msmsScans.txt", sep = ''))
  msscans <- msscans[, c('Raw file', 'Total ion current', 'Base peak intensity', 'Identified', 'MS/MS IDs', 'Sequence', 'Proteins', 'Experiment')]
  msscans <- msscans[msscans$Identified == '+' ,]
  msscans <- msscans[grep(paste(isoform1, isoform2, sep = '|'), msscans$Proteins, value = F), ]
  
  evidence <- fread(paste(path.to.mq, 'evidence.txt', sep = ''))
  peptides.both <- as.character(unique(evidence[grep(paste(isoform1, ".*", isoform2, "|", 
                                                           isoform2, ".*", isoform1, sep = ''),
                                                     evidence$Proteins), c('Sequence')]))
  peptides.isoform1 <- as.character(intersect(unique(evidence[grep(isoform1, evidence$Proteins), c('Sequence')]), 
                                 unique(evidence[grep(isoform2, evidence$Proteins, invert = T), c('Sequence')])))
  peptides.isoform2 <- as.character(intersect(unique(evidence[grep(isoform2, evidence$Proteins), c('Sequence')]), 
                                              unique(evidence[grep(isoform1, evidence$Proteins, invert = T), c('Sequence')])))
  msscans$protein <- 'both'  
  msscans[msscans$Sequence %in% peptides.isoform1,]$protein <- 
    as.character(intersect(unique(evidence[grep(isoform1, evidence$Proteins), c('Gene names')]), 
                           unique(evidence[grep(isoform2, evidence$Proteins, invert = T), c('Gene names')])))
  msscans[msscans$Sequence %in% peptides.isoform2,]$protein <- 
    as.character(intersect(unique(evidence[grep(isoform2, evidence$Proteins), c('Gene names')]), 
                           unique(evidence[grep(isoform1, evidence$Proteins, invert = T), c('Gene names')])))
  evidence <- NULL
  msscans <- msscans[, c('Experiment', 'protein', 'Total ion current')]
  msscans <- aggregate(msscans, by = list(msscans$Experiment, msscans$protein), FUN = mean)
  msscans <- msscans[, c("Group.1", "Group.2", "Total ion current")]  
  colnames(msscans) <- c('experiment', 'protein', 'Total ion current')
  if(old.data){
    msscans$experiment <- gsub('^NCBP', 'Rockefeller_NCBP', msscans$experiment)
    msscans$experiment <- gsub('^tk', 'Columbia', msscans$experiment)
    msscans$target <- 'Rockefeller NCBP3'
    msscans[grep('NCBP1', msscans$experiment), ]$target <- 'Rockefeller NCBP1'
    msscans[grep('NCBP2', msscans$experiment), ]$target <- 'Rockefeller NCBP2'
    msscans[grep('Columbia_NCBP3', msscans$experiment), ]$target <- 'Columbia NCBP3'
    msscans$target <- factor(msscans$target, levels = c('Rockefeller NCBP1', 'Rockefeller NCBP2', 'Rockefeller NCBP3',
                                                        'Columbia NCBP3'))
    msscans$experiment <- factor(msscans$experiment, levels = c(grep('Rockefeller_NCBP1', unique(msscans$experiment),
                                                                     value = T),
                                                                grep('Rockefeller_NCBP2', unique(msscans$experiment),
                                                                     value = T),
                                                                grep('Rockefeller_NCBP3', unique(msscans$experiment),
                                                                     value = T),
                                                                grep('Columbia_NCBP3', unique(msscans$experiment),
                                                                     value = T)))
  } else{
    pretty.names <- fread('data/pretty_experiments_names.csv')
    pretty.names$target <- str_split_fixed(pretty.names$experiment, '_', 2)[, 2]
    msscans$target <- str_split_fixed(msscans$experiment, '_', 2)[, 1]
    msscans$location <- 'Rockefeller'
    msscans[msscans$target == 'NCBP3', ]$location <- 'Columbia'
    pretty.names$location <- str_split_fixed(pretty.names$experiment, '_', 2)[, 1]
    
    msscans$experiment2 <-  unname(unlist(apply(msscans, 1, function(x) unique(pretty.names[pretty.names$target == x['target'] & 
                                                                                              pretty.names$location == x['location'], ]$pretty_name))))
    msscans$buffer <- str_split_fixed(msscans$experiment, '_', 2)[, 2]
    msscans$experiment <- paste(msscans$experiment2, msscans$buffer, sep = '_')
  }
  pdf(paste("out/QC/", MakeDate(), "_isoforms_merging.pdf", sep = ''), width = 13, height = 10)
  print(ggplot(msscans, aes(x = experiment, y = log(`Total ion current`), col = protein, group = protein)) + 
    geom_point() + 
    geom_line() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_grid(rows = vars(target), scales = "free_x"))
  dev.off()
  return(0)
  
}
  





