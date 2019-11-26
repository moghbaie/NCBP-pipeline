source('src/MaxQuant.R')
source('lib/constants.R')
source('src/Preprocessing/preprocessing.R')
source('src/Imputation/impute.R')

library(stringr)
library(reshape2)
library(gridExtra)
library(grid)
require(stringr)

PlotImputeGroup <- function(intensities.group, experiments.case, 
                            experiments.control, do.plot = T, path., metadata.experiments.) {
  
  group.name.old <- gsub('_', '.', experiments.case[1], fixed = T)
  group.name.old <- substr(group.name.old, 1, nchar(group.name.old) - 1)
  group.name <- 
    paste(metadata.experiments.[metadata.experiments.$experiment.name == experiments.case[1], 
                                c('pretty_name', 'buffer'), with = F], collapse = ', condition ')
  print(group.name)
  print('case')
  intensities.group.case <- intensities.group[, c(experiments.case, 'Protein IDs'), with = F]
  plots.case <- list()
  files.list <- list.files(path = path., pattern = 'imputed', full.names = T)
  
  for (method in names(kOldImputationMethods)) {
    # print(method)
    
    files.found <- grep(paste(kOldImputationMethods[[method]]$sep.method,
                              kOldImputationMethods[[method]]$method.mar,
                              kOldImputationMethods[[method]]$method.mnar, sep = "_"),
                        files.list, value = T)
    
    intensities.imputed <- readRDS(files.found[length(files.found)])
    intensities.imputed <- intensities.imputed[, c(experiments.case, 'Protein IDs'), with = F]
    intensities.imputed <- intensities.imputed[complete.cases(intensities.imputed), ] 
    intensities.imputed <- intensities.imputed[order(intensities.imputed$`Protein IDs`), ]
    
    merged <- merge(intensities.imputed, intensities.group.case, by = "Protein IDs", 
                    suffixes = c(".imp",".old"))
    merged$n.zeros <- rowSums(merged == 0)
    merged <- merged[order(merged$n.zeros, decreasing = T), ]
    
    corr.m <- round(cor(merged[, 
                               paste(colnames(intensities.imputed)[colnames(intensities.imputed) 
                                                                   %nin% "Protein IDs"], '.imp',
                                     sep = ''),
                               with = F]), 2)
    
    rownames(corr.m) <- gsub('.imp', '', gsub(group.name.old, '', rownames(corr.m)))
    colnames(corr.m) <- gsub('.imp', '', gsub(group.name.old, '', colnames(corr.m)))
    melted.cormat <- reshape2::melt(corr.m)
    
    merged$take <- F
    merged <- merged[, c(2, 3, ncol(intensities.group.case) + 1, ncol(intensities.group.case) + 2,
                         (ncol(merged) - 1), ncol(merged)), with = F]
    colnames(merged) <- c('rep.a.imp', 'rep.b.imp', 'rep.a.old', 'rep.b.old', 'n.zeros', 'take')
    merged[merged$n.zeros == (ncol(intensities.group.case) - 1),]$take <- T
    
    if (!empty(merged[merged$n.zeros == 0,])) {
      merged[merged$n.zeros == 0,]$take <- T
    }
    
    merged$take <- T
    merged <- merged[merged$take == T,]
    total <- data.frame('rep.a' = c(merged$rep.a.old, merged$rep.a.imp), 
                        'rep.b' = c(merged$rep.b.old, merged$rep.b.imp))
    total$impute <- c(rep('no', nrow(merged)), rep('yes', nrow(merged)))
    total$n.zeros <- c(rep(4, nrow(merged[merged$n.zeros == 4, ])), 
                       rep(3, nrow(merged[merged$n.zeros == 3, ])), 
                       rep(2, nrow(merged[merged$n.zeros == 2, ])),
                       rep(1, nrow(merged[merged$n.zeros == 1, ])),
                       rep(0, nrow(merged[merged$n.zeros == 0, ])))
    total$n.zeros <- as.factor(total$n.zeros)
    total$impute <- as.factor(total$impute)
    
    p <- ggplot(total, aes(x = rep.a, y = rep.b, color = n.zeros, shape = impute)) + 
      geom_point(size = 3) + 
      geom_abline(slope = 1, intercept = 0) + 
      theme(plot.title = element_text(size = 8, face = "bold")) + 
      theme(legend.spacing.y = unit(0.1, "cm")) +
      ggtitle(label = paste(gsub('method', 'step 3.', str_split(method, pattern = "\\.")[[1]][2]), 
                            paste0('step 4.', str_split(method, pattern = "\\.")[[1]][1]), 
                            sep = ', ')) +
      labs(color = "# zero replicates") +
      scale_shape_discrete(name = "", labels = c("before imputation", "after imputation"))
    
    plots.case[[method]] <- p
    
    if (!is.na(str_split(method, pattern = "\\.")[[1]][2]) & 
        str_split(method, pattern = "\\.")[[1]][2] == '3') {
      plots.case[[paste(method, 'corr', sep = '.')]] <- 
        ggplot(data = melted.cormat, aes(x = Var1, y = Var2, fill = value)) + 
        geom_tile() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = ' ', y = ' ')
    }
  }
  
  intensities.group.control <- intensities.group[, c(experiments.control, 'Protein IDs'), with = F]
  plots.control <- list()
  
  print('control')
  for (method in names(kOldImputationMethods)) {
    # print(method)
      
    files.found <- grep(paste(kOldImputationMethods[[method]]$sep.method,
                              kOldImputationMethods[[method]]$method.mar,
                              kOldImputationMethods[[method]]$method.mnar, sep = "_"),
                        files.list, value = T)
    
    intensities.imputed <- readRDS(files.found[length(files.found)])
      intensities.imputed <- intensities.imputed[, c(experiments.control, 'Protein IDs'), with = F]
      intensities.imputed <- intensities.imputed[complete.cases(intensities.imputed), ] 
      intensities.imputed <- intensities.imputed[order(intensities.imputed$`Protein IDs`), ]    
      merged <- merge(intensities.imputed, intensities.group.control, by = "Protein IDs", 
                      suffixes = c(".imp",".old"))
      merged$n.zeros <- rowSums(merged == 0)
      merged <- merged[order(merged$n.zeros, decreasing = T), ]
      corr.m <- 
        round(cor(merged[, paste(colnames(intensities.imputed)[colnames(intensities.imputed) 
                                                          %nin% "Protein IDs"], '.imp', sep = ''),
                         with = F]), 2)
      
      rownames(corr.m) <- gsub(paste(unlist(strsplit(group.name.old, '.', fixed = T))[1], 
                                     unlist(strsplit(group.name.old, '.', fixed = T))[2], "C", 
                                     unlist(strsplit(group.name.old, '.', fixed = T))[3], sep = '_'),
                               '', rownames(corr.m))
      colnames(corr.m) <- gsub(paste(unlist(strsplit(group.name.old, '.', fixed = T))[1], 
                                     unlist(strsplit(group.name.old, '.', fixed = T))[2], "C", 
                                     unlist(strsplit(group.name.old, '.', fixed = T))[3], sep = '_'),
                               '', colnames(corr.m))
      
      rownames(corr.m) <- gsub('.imp', '', gsub(group.name.old, '', rownames(corr.m)))
      colnames(corr.m) <- gsub('.imp', '', gsub(group.name.old, '', colnames(corr.m)))
      melted.cormat <- reshape2::melt(corr.m)
      merged$take <- F
      merged <- merged[, c(2, 3, ncol(intensities.group.control) + 1, 
                           ncol(intensities.group.control) + 2,
                           (ncol(merged) - 1), ncol(merged)), with = F]
      
      colnames(merged) <- c('rep.a.imp', 'rep.b.imp', 'rep.a.old', 'rep.b.old', 'n.zeros', 'take')
      
      merged[merged$n.zeros == (ncol(intensities.group.control) - 1),]$take <- T
      
      if (!empty(merged[merged$n.zeros == 0,])) {
        merged[merged$n.zeros == 0,]$take <- T
      }
      
      merged$take <- T
      merged <- merged[merged$take == T,]
      total <- data.frame('rep.a' = c(merged$rep.a.old, merged$rep.a.imp), 
                          'rep.b' = c(merged$rep.b.old, merged$rep.b.imp))
      
      total$impute <- c(rep('no', nrow(merged)), rep('yes', nrow(merged)))
      total$n.zeros <- c(rep(4, nrow(merged[merged$n.zeros == 4, ])), 
                         rep(3, nrow(merged[merged$n.zeros == 3, ])), 
                         rep(2, nrow(merged[merged$n.zeros == 2, ])),
                         rep(1, nrow(merged[merged$n.zeros == 1, ])),
                         rep(0, nrow(merged[merged$n.zeros == 0, ])))
      
      total$n.zeros <- as.factor(total$n.zeros)
      total$impute <- as.factor(total$impute)
      
      p <- ggplot(total, aes(x = rep.a, y = rep.b, color = n.zeros, shape = impute)) + 
        geom_point(size = 3) + 
        geom_abline(slope = 1, intercept = 0) + 
        theme(plot.title = element_text(size = 8, face = "bold")) + 
        theme(legend.spacing.y = unit(0.1, "cm")) +
        ggtitle(label = paste(gsub('method', 'step 3.', str_split(method, pattern = "\\.")[[1]][2]), 
                                    paste0('step 4.', str_split(method, pattern = "\\.")[[1]][1]), 
                              sep = ', ')) +
        labs(color = "# zero replicates") +
        scale_shape_discrete(name = "", labels = c("before imputation", "after imputation"))
      
      plots.control[[method]] <- p
      
      if (!is.na(str_split(method, pattern = "\\.")[[1]][2]) & 
          str_split(method, pattern = "\\.")[[1]][2] == '3') {
        plots.control[[paste(method, 'corr', sep = '.')]] <- 
          ggplot(data = melted.cormat, aes(x = Var1, y = Var2, fill = value)) + 
          geom_tile() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(x = '', y = '')
      }
  }
  
  if (do.plot) 
    pdf(paste('out/', group.name, 'imputation.replicates.correlation.pdf', sep = ''), height = 15, 
        width = 10)
  
  do.call("grid.arrange", c(plots.case, ncol = 4, 
                            top = paste(group.name, '\n cases', sep = '')))
  do.call("grid.arrange", c(plots.control, ncol = 4, 
                            top = paste(group.name, '\n controls', sep = '')))
  
  if (do.plot) 
    dev.off()
} 



PlotImputeAll <- function(max.quant.obj, path.) {
  date <- MakeDate()
  
  pdf(paste(path., date, '_replicates.correlation.pdf', sep = ''), height = 20, width = 15)
  
  intensities.before.impute <- max.quant.obj$GetIntensityTable()
  metadata.experiments <- max.quant.obj$GetMetadataExperiments()
  metadata.experiments$group <- paste(metadata.experiments$location, 
                                          metadata.experiments$target, 
                                          metadata.experiments$buffer, sep = '_')
  metadata.experiments$buffer <- as.numeric(metadata.experiments$buffer)
  setorderv(metadata.experiments, c("location", "target", "buffer"), order = c(-1, 1, 1))
  
  for (exp.group in unique(metadata.experiments$group)) {
    
    experiments <- metadata.experiments[metadata.experiments$group == exp.group,]$experiment.name
    intensities.group <- intensities.before.impute[, c(experiments, 'Protein IDs'), with = F]
    intensities.group <- intensities.group[rowSums(intensities.group[, experiments, with = F]) != 0,] 
    intensities.group <- intensities.group[order(intensities.group$`Protein IDs`),]
    
    experiments.case = 
      metadata.experiments[metadata.experiments$group == exp.group &
                             metadata.experiments$is.case == T,]$experiment.name
    experiments.control = 
      metadata.experiments[metadata.experiments$group == exp.group &
                             metadata.experiments$is.case == F,]$experiment.name
    
    PlotImputeGroup(intensities.group, experiments.case, 
                    experiments.control, do.plot = F,
                    path. = path., 
                    metadata.experiments. = metadata.experiments)
  }
  
  dev.off()
  return(0)
}

CheckReplicatesCorrelationAfterImputation <- function(path = "./data/txt_noGA/proteinGroups.txt") {
  
  for (method in names(kOldImputationMethods)) {
      print(method)
    
      max.quant <- BuildMaxQuant(path, pattern = pattern,
                                 groups.list = kGroupsList, subgroup.list = kSubgroupList)
      
      max.quant$RemoveContaminantsReversed()$
        SumMAGOH(save.result.table = F)$
        LogTricky()
      
      max.quant$ImputeAll(
        id.column = "Protein IDs",
        sep.method = kOldImputationMethods[[method]]$sep.method,
        method.mar = kOldImputationMethods[[method]]$method.mar,
        method.mnar = kOldImputationMethods[[method]]$method.mnar,
        sigma.const = 4,
        save.result.table = T,
        path. = 'out/imputation/'
      )
  }
  
  max.quant <- BuildMaxQuant(path, pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList)
  max.quant$MakePrettyExperimentNames()
  max.quant$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky()
  PlotImputeAll(max.quant, path. = 'out/imputation/')
}


CompareWithPerseus <- function(path.to.perseus.output) {
  # manual
  max.quant <- BuildMaxQuant(paste0(path.to.perseus.output,"/proteinGroups.txt"), 
                             pattern = "LFQ intensity",
                             groups.list = kGroupsList, subgroup.list = kSubgroupList)
  metadata.experiments <- max.quant$GetMetadataExperiments()
  metadata.experiments$replicate <- metadata.experiments$buffer
  metadata.experiments[grep('Control', metadata.experiments$experiment.name), ]$buffer <- 
    str_split_fixed(metadata.experiments[grep('Control', 
                                              metadata.experiments$experiment.name), 
                                         ]$experiment.name, '_', 3)[, 2]  
  metadata.experiments[grep('Control', metadata.experiments$experiment.name), ]$is.case <- FALSE
  metadata.experiments[grep('Control', metadata.experiments$experiment.name), ]$target <- 'NCBP3'
  metadata.experiments[grep('Control', metadata.experiments$experiment.name), ]$location <- 'Columbia'
  metadata.experiments[grep('Control', 
                            metadata.experiments$experiment.name, invert = T), ]$buffer <- 
    str_split_fixed(metadata.experiments[grep('Control', metadata.experiments$experiment.name, 
                                              invert = T), ]$experiment.name, '_', 4)[, 3]  
  max.quant$SetMetadataExperiments(metadata.experiments, change.string = 'Fixed metadata')
  
  max.quant$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky()
  
  max.quant$ImputeAll(id.column = "Protein IDs", sep.method = "all.mnar",
                      method.mar = "",
                      method.mnar = 'perseus', 
                      sigma.const = 4, save.result.table = F)
  
  
  # perseus
  date <- paste(str_split_fixed(Sys.time(), ' ', 3)[, 1], str_split_fixed(Sys.time(), ' ', 3)[, 2], 
                sep = '_')
  pdf(paste('out/imputation/', date, '_compare.perseus.our.perseus.pdf', sep = ''))
  
  for (buffer in unique(max.quant$GetMetadataExperiments()$buffer)) {
    print(buffer)
    our.perseus <-
      max.quant$GetIntensityTable()[, c('Protein IDs', 
                                        grep(buffer, colnames(max.quant$GetIntensityTable()), 
                                                          value = T)), with = F]
    perseus <- fread(paste(path.to.perseus.output, '/', buffer, '.txt', sep = ''), header = TRUE, 
                     sep = "\t")
    perseus <- perseus[-c(1, 2), ]
    lfq.columns <- grep(x = colnames(perseus), pattern = "^LFQ intensity", value = T)
    protein.groups.columns <- c("Protein IDs", lfq.columns)
    perseus <- perseus[, protein.groups.columns, with = F]
    setnames(perseus, lfq.columns, gsub(".* (.*)", "\\1", lfq.columns))
    our.perseus <- max.quant$GetIntensityTable()
    our.perseus <- our.perseus[, c(grep(buffer, colnames(our.perseus), value = T), 'Protein IDs'), with = F]
    colnames(our.perseus) <- gsub('Columbia_', '', colnames(our.perseus))
    merged <- merge(perseus, our.perseus, by = 'Protein IDs', suffixes = c('.perseus','.our.perseus'))
    
    for (replica in unique(max.quant$GetMetadataExperiments()$replicate)){
      for (is.case in c(T, F)){
        if (is.case) {
          MIN <- round(min(as.numeric(unlist(c(merged[, paste0('NCBP3_', buffer, '_', 
                                                               replica, '.perseus'), with = F], 
                                               merged[, paste0('NCBP3_', buffer, '_', 
                                                               replica, '.our.perseus'), with = F])))))
          MAX <- round(max(as.numeric(unlist(c(merged[, paste0('NCBP3_', buffer, '_', 
                                                               replica, '.perseus'), with = F], 
                                               merged[, paste0('NCBP3_', buffer, '_', 
                                                               replica, '.our.perseus'), with = F])))))
          plot(x = unlist(merged[, paste0('NCBP3_', buffer, '_', replica, '.perseus'), with = F]),
             y = unlist(merged[, paste0('NCBP3_', buffer, '_', replica, '.our.perseus'), with = F]), 
             main = paste0('NCBP3-LAP, test tubes, condition ', buffer, ', replica ', replica, ', is_case ', is.case),
             xlab = 'perseus', ylab = 'manual implementation', xlim = c(MIN, MAX), ylim = c(MIN, MAX), axes = F)
          axis(side = 1, at = seq(MIN, MAX, 5))
          axis(side = 2, at = seq(MIN, MAX, 5))
        } else { 
          MIN <- round(min(as.numeric(unlist(c(merged[, paste0('Control_', buffer, '_', replica, '.perseus'), with = F], 
                                               merged[, paste0('Control_', buffer, '_', replica, '.our.perseus'), with = F])))))
          MAX <- round(max(as.numeric(unlist(c(merged[, paste0('Control_', buffer, '_', replica, '.perseus'), with = F], 
                                               merged[, paste0('Control_', buffer, '_', replica, '.our.perseus'), with = F])))))
          plot(x = unlist(merged[, paste0('Control_', buffer, '_', replica, '.perseus'), with = F]),
               y = unlist(merged[, paste0('Control_', buffer, '_', replica, '.our.perseus'), with = F]),
               main = paste0('NCBP3-LAP, test tubes, condition ', buffer, ', replica ', replica, ', is_case ', is.case),
               xlab = 'perseus', ylab = 'manual implementation', xlim = c(MIN, MAX), ylim = c(MIN, MAX), axes = F)
          axis(side = 1, at = seq(MIN, MAX, 5))
          axis(side = 2, at = seq(MIN, MAX, 5))
        } 
      }
    }
  }
  dev.off()  
}

