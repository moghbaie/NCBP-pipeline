library(data.table)
library(stringr)
library(bit64)
library(scales)
library(parallel)
library(pbmcapply)
source("lib/constants.R")
source("src/Imputation/set.null.and.impute.R")


CountRMSEsOnce <- function(max.quant.obj, sigma.const = 4, method.allnulls = "method.1",
                           seed = seed) {
  results <- list()
  rmse.methods <- list()
  square.errors.methods <- c()
  
  nrmses <- c()
  maxs <- c()
  mins <- c()
  intensities <- max.quant.obj$GetIntensityTable()
  
  for(colname in colnames(intensities)[which(colnames(intensities) != 'Protein IDs')]){
    maxs <- c(maxs, max(intensities[which(intensities[, colname, with = F] != 0), colname, with = F]))
    mins <- c(mins, min(intensities[which(intensities[, colname, with = F] != 0), colname, with = F]))
  }
  
  names(mins) <- colnames(intensities)[which(colnames(intensities) != 'Protein IDs')]
  names(maxs) <- colnames(intensities)[which(colnames(intensities) != 'Protein IDs')]
  mins.type <- list()
  maxs.type <- list()
  
  metadata <- max.quant.obj$GetMetadataExperiments()
  
  if (is.null(metadata$group)){
    metadata$group <- paste(metadata$location, metadata$target, metadata$buffer, sep = '_')
  } 
  
  groups <- unique(metadata$group)
  
  for (exp in groups) {
    for (is.case. in c(T, F)) {
      name <- metadata[metadata$group == exp & metadata$is.case == is.case.,]$experiment.name[1]
      mins.type[[substr(name, 1, nchar(name) - 1)]] <- min(mins[metadata[metadata$group == exp & 
                                                                           metadata$is.case == is.case., ]$experiment.name])
      maxs.type[[substr(name, 1, nchar(name) - 1)]] <- max(maxs[metadata[metadata$group == exp & 
                                                                           metadata$is.case == is.case., ]$experiment.name])
    }
  }
  
  nrmses <- setNames(data.frame(matrix(ncol = length(groups), nrow = length(methods.notnulls) + 1)), groups)
  rownames(nrmses) <- c(methods.notnulls, 'random')
  
  for (method.not.nulls in methods.notnulls) {
    square.errors.method <- c()
    square.errors.random <- c()
    print(method.not.nulls)
    
    if (method.not.nulls == "perseus") {
      old.method.name = "perseus"
    } else {
      old.method.name <- paste0(method.not.nulls, ".", str_remove(method.allnulls, ".*\\.")) 
    }

    for (exp in groups) {
      for (is.case. in c(T, F)) {
        cols.to.use <-
          metadata[metadata$group == exp & metadata$is.case == is.case.,]$experiment.name
        list.type <-
          MakeRandomAndImputedMtx(
            intensities = max.quant.obj$GetIntensityTable(),
            cur.exp.cols = cols.to.use,
            sep.method = kOldImputationMethods[[old.method.name]]$sep.method,
            method.mar = kOldImputationMethods[[old.method.name]]$method.mar,
            method.mnar = kOldImputationMethods[[old.method.name]]$method.mnar,
            sigma.const = sigma.const,
            random.method = 'rnorm',
            seed = seed)
        
        to.remove <- which(list.type$rows.deleted %in% list.type$rows.all.nulls)
        rows.to.take <- list.type$rows.deleted[-to.remove]
        cols.to.take <- list.type$cols.deleted[-to.remove]
        
        
        square.errors.type <- mapply(function(row, column) (list.type$original[row, column] -
                 list.type$imputed[row, column]) ^ 2, rows.to.take, cols.to.take)
        square.errors.random <- c(square.errors.random, mapply( function(row, column)
                  (list.type$original[row, column] - list.type$random[row, column]) ^ 2, 
                  rows.to.take, cols.to.take))
        square.errors.method <- c(square.errors.method, square.errors.type)
        
        name <- metadata[metadata$group == exp & metadata$is.case == is.case.,]$experiment.name[1]
        nrmse <- sqrt(mean(square.errors.type, na.rm = T)) / (maxs.type[[substr(name, 1, (nchar(name) - 1))]] - 
                                                                mins.type[[substr(name, 1, (nchar(name) - 1))]])
        nrmse.random <- sqrt(mean(square.errors.random, na.rm = T)) / (maxs.type[[substr(name, 1, (nchar(name) - 1))]] - 
                                                                         mins.type[[substr(name, 1, (nchar(name) - 1))]])
        nrmses[method.not.nulls, substr(name, 1, (nchar(name) - 1))] <- nrmse
        nrmses['random', substr(name, 1, (nchar(name) - 1))] <- nrmse.random
      }
    }
    rmse.methods[[paste(method.not.nulls, method.allnulls, sep = '.')]] <- sqrt(mean(square.errors.method, na.rm = T))
    square.errors.methods[[paste(method.not.nulls, method.allnulls, sep = '.')]] <- square.errors.method
  }
  
  rmse.random <- sqrt(mean(square.errors.random, na.rm = T))
  rmses <- as.data.frame(unlist(rmse.methods))
  colnames(rmses) <- 'rmse'
  rmses <- rbind(rmses, 'random' = c(rmse.random))
  results <- list("rmses" = rmses, 'errors' = square.errors.methods, 'nrmses' = nrmses)
  return(results)
}


CountRMSEConfIntervals <- function(rmses){
  saveRDS(rmses, paste0('out/imputation/', MakeDate(), 'rmses.rds'))
  stat <- rbind('means' = colMeans(rmses), 'sds' = sapply(data.frame(rmses), sd)) 
  stat <- as.data.frame(t(stat))
  stat <- stat[order(stat$means), ]
  conf.intervals <- data.frame('method' = rownames(stat), 
                               'rmse' = paste(round(stat$means, 2), '±', 
                                              round(qt(0.975, df = nrow(rmses) - 1) *
                                                      sqrt(stat$sds / nrow(rmses)), 2), sep = ''))
  conf.intervals$method <- as.character(conf.intervals$method)
  conf.intervals$method <- gsub('.method.1', '', conf.intervals$method)
  conf.intervals$method <- gsub('method', '', conf.intervals$method)
  
  colnames(conf.intervals) <- c('step4', 'RMSE')
  
  pdf(paste0('out/imputation/', MakeDate(), '.conf.interval.rmse.with.pdf'), width = 3, height = 5)
  grid.table(conf.intervals, rows = NULL)
  dev.off()
}

PlotNRMSE <- function(nrmses, metadata.experiments){
  saveRDS(nrmses, paste0('out/imputation/', MakeDate(), 'nrmses.rds'))
  x <- lapply(nrmses, function(dt){
    dt$method <- row.names(dt)
    melted <- melt(dt, id.vars = 'method', variable.name = 'experiment')
    melted
  })
  nrmses.dt <- as.data.table(Reduce(function(...) merge(..., by = c('method', 'experiment')), x))
  nrmses.dt$mean. <- rowMeans(nrmses.dt[, colnames(nrmses.dt) %nin% c('method', 'experiment'), with = F])
  nrmses.dt$sds <- unlist(apply(nrmses.dt[, colnames(nrmses.dt) %nin% c('method', 'experiment', 'mean.'), with = F], 1, sd))
  nrmses.dt$error <- unlist(apply(nrmses.dt[, colnames(nrmses.dt) %nin% c('method', 'experiment', 'mean.'), with = F], 1, function(row) {
    round(qt(0.975, df = length(row) - 2) * sqrt(row['sds'] / length(row)), 2)
  }))
  colours <- unname(kColorList[2:(length(unique(nrmses.dt$method)) + 1)])
  nrmses.dt$method <- gsub('method', 'step4.', nrmses.dt$method)
  
  nrmses.dt$experiment <- unlist(lapply(nrmses.dt$experiment, 
                                        function(experiment.) {
                                          x <- paste(metadata.experiments[metadata.experiments$experiment.name == paste0(experiment., 'a'), 
                                                                          c('pretty_name', 'buffer'), with = F], collapse = ', condition ')
                                          if(metadata.experiments[metadata.experiments$experiment.name == paste0(experiment., 'a'), 
                                                                  c('is.case'), with = F]$is.case == F) {x <- paste0(x, ', Control')}
                                          x
                                        }))
  nrmses.dt <- nrmses.dt[, c("method", "experiment", "mean.", "sds", "error"), with = F]
  pdf(paste('out/imputation/', MakeDate(), '_nrmse.pdf', sep = ''), width = 13, height = 8)
  pd <- position_dodge(0.1)
  p <- ggplot(nrmses.dt, aes(x = experiment, y = mean., colour = as.factor(method), group = method)) + 
    geom_point() + 
    geom_line() + 
    geom_errorbar(aes(ymin = mean. - error, ymax = mean. + error), 
                  width = .1, position = pd) + 
    ylab(label = 'NRMSE') + 
    scale_y_log10() +
    ggtitle(label = 'Normalized Root-mean-square error for methods in various experiments') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    scale_color_manual(values = unname(colours)) + 
    theme(legend.text = element_text(size = 15)) + 
    guides(colour = guide_legend(title = ""))
  print(p)
  dev.off()  
}

CountRMSES <- function(max.quant.obj,
                       n.threads = 7,
                       n.iterations = 20){
  global.rmses <- pbmclapply(c(1:n.iterations), function(seed.){
    results <- CountRMSEsOnce(max.quant.obj)
    rmses <- t(results$rmses)
    rmses <- rmses[, order(colnames(rmses))]
    nrmses <- results$nrmses 
    list(rmses, nrmses)
  }, mc.cores = n.threads, mc.silent = F)
  
  rmses <- do.call('rbind', lapply(global.rmses, function(result) result[[1]])) 
  nrmses <- lapply(global.rmses, function(result) result[[2]])
  
  CountRMSEConfIntervals(rmses)
  PlotNRMSE(nrmses, max.quant.obj$GetMetadataExperiments())
  
}



