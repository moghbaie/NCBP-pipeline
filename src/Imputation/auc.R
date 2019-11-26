library(data.table)
library(stringr)
library(bit64)
library(gplots)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(parallel)
library(pbmcapply)

source("lib/constants.R")
source("src/Imputation/set.null.and.impute.R")
source("lib/helper_functions.R")

MakeRocCurve <- function(deltas.random, deltas.imputed, title = '', subtitle = '', 
                           xlab = 'delta', ylab="y", reverse = T, col = "red", 
                           add = F, plot.hist = F) {
  
  hist.random <- hist(deltas.random, col = scales::alpha('red',.5), border=F, xlim = c(-1, 1),
                      ylim = c(0, 40), breaks = seq(from = - 2, to = 2, by = 0.02),
                      plot = plot.hist, xlab = xlab, main = title)
  hist.imputed <- hist(deltas.imputed, add = T, plot = plot.hist,  border=F, seq(from = -2, 
                                                                               to = 2, by = 0.02), col = scales::alpha('skyblue', .6))
  if (plot.hist) {
    legend("topleft", bty = "n", fill = c("skyblue", "red"), legend = c("imputed", "random"))
  }
  
  # making ROC curve
  if (reverse) {
    roc.y <- cumsum(rev(hist.random$counts / sum(hist.random$counts) * 100))
    roc.x <- cumsum(rev(hist.imputed$counts / sum(hist.imputed$counts) * 100))
  } else {
    roc.y <- cumsum(hist.random$counts / sum(hist.random$counts) * 100)
    roc.x <- cumsum(hist.imputed$counts / sum(hist.imputed$counts) * 100)
  }
  
  dt.roc = data.table(x = roc.x, y = roc.y)
  auc = 0
  for (i in c(1:(length(roc.x) - 1))) {
    if ((roc.x[i + 1] - roc.x[i]) > 0) {
      auc = auc + (roc.y[i] + roc.y[i + 1]) / 2 * (roc.x[i + 1] - roc.x[i])
    }
  }
  
  if (!plot.hist) {
    if(!add){
      plot(dt.roc, type = "l", frame = F, col = col, xlab = xlab, ylab = ylab, 
           main = title, lwd = 2)
      abline(a = 0, b = 1, col = "blue", lty = 2)
    } else {
      lines(dt.roc, col = col, xlab = xlab, ylab = ylab, lwd = 2)
    }
  }
  
  auc / 10000
}

MakeRocCurveDensity <- function(cur.exp.cols, method.not.nulls, method.allnulls, 
                                   is_log_tranformed = T, sigma_const = 3, title) {
  
  mtx.list <- MakeRandomAndImputedMtx(cur.exp.cols, method.not.nulls, method.allnulls, 
                                          is_log_tranformed = T, sigma_const = 3) 
  
  deltas.random <- as.vector((mtx.list$original - mtx.list$random)*2/(mtx.list$original +
                                                                        mtx.list$random))
  deltas.imputed <- as.vector((mtx.list$original - mtx.list$imputed)*2/(mtx.list$original + 
                                                                          mtx.list$imputed))
  
  deltas.random <- deltas.random[deltas.random != 0]
  deltas.imputed <- deltas.imputed[deltas.imputed != 0]
  
  MakeRocCurve(deltas.random = deltas.random, deltas.imputed = deltas.imputed,
                 title = title, subtitle = paste0(method.not.nulls, ", ", method.allnulls))
}


CountCorrRandomImputed <- function(intensities, sep.method, method.mar, method.mnar, 
                                   use_only_not_nulls, sigma_const, metadata,
                                   compare.with.mnar.1.nonzero = F, seed) {
  if (!missing(seed))
    set.seed(seed)
  
  all.corr.rand <- c()
  all.corr.imputed <- c()
  
  metadata$group <- paste(metadata$location, metadata$target, metadata$buffer, sep = '_')
  
  #controls:
  print('controls')
  
  for (exp.group in unique(metadata$group)) {
    print(exp.group)
    cur.exp.cols <- metadata[metadata$group == exp.group & metadata$is.case == F, ]$experiment.name
    mtx.list <- MakeRandomAndImputedMtx(intensities, cur.exp.cols, 
                                        sep.method, method.mar, method.mnar, 
                                        sigma_const,
                                        compare.with.mnar.1.nonzero = compare.with.mnar.1.nonzero,
                                        seed = seed) 
    
    for (i in c(1:length(cur.exp.cols))) {
      if (use_only_not_nulls) {
        to_remove <- which(mtx.list$rows.deleted %in% mtx.list$rows.all.nulls)
        rows.to_take <- mtx.list$rows.deleted[-to_remove]
        cols.to_take <- mtx.list$cols.deleted[-to_remove]
        rows.were.not.null <- rows.to_take[which(cols.to_take == i)]  
        
      } else {
        rows.were.not.null <- mtx.list$rows.deleted[which(mtx.list$cols.deleted == i)]  
      }
      
      
      all.corr.rand <- c(all.corr.rand, cor(mtx.list$original[rows.were.not.null, i],
                                            mtx.list$random[rows.were.not.null, i]))
      all.corr.imputed <- c(all.corr.imputed, cor(mtx.list$original[rows.were.not.null, i],
                                                  mtx.list$imputed[rows.were.not.null, i]))
    }    
  }

  
  #cases:
  print('cases')
  
  for (exp.group in unique(metadata$group)) {
    print(exp.group)
    cur.exp.cols <- metadata[metadata$group == exp.group & metadata$is.case == T, ]$experiment.name
    mtx.list <- MakeRandomAndImputedMtx(intensities, cur.exp.cols, 
                                        sep.method, method.mar, method.mnar, 
                                        sigma_const,
                                        compare.with.mnar.1.nonzero = compare.with.mnar.1.nonzero,
                                        seed = seed) 
    
    for (i in c(1:length(cur.exp.cols))) {
      if (use_only_not_nulls) {
        to_remove <- which(mtx.list$rows.deleted %in% mtx.list$rows.all.nulls)
        rows.to_take <- mtx.list$rows.deleted[-to_remove]
        cols.to_take <- mtx.list$cols.deleted[-to_remove]
        rows.were.not.null <- rows.to_take[which(cols.to_take == i)]  
        
      } else {
        rows.were.not.null <- mtx.list$rows.deleted[which(mtx.list$cols.deleted == i)]  
      }
      
      all.corr.rand <- c(all.corr.rand, cor(mtx.list$original[rows.were.not.null, i],
                                            mtx.list$random[rows.were.not.null, i]))
      all.corr.imputed <- c(all.corr.imputed, cor(mtx.list$original[rows.were.not.null, i],
                                                  mtx.list$imputed[rows.were.not.null, i]))
    }
  }
  
  list(corr.rand = all.corr.rand, corr.imputed = all.corr.imputed)
}

PlotROCcurves <- function(max.quant.obj, method.allnulls = "method.1") {
  intensities <- max.quant.obj$GetIntensityTable()
  metadata <- max.quant.obj$GetMetadataExperiments()
  pdf(paste0('out/imputation/', MakeDate(), "_roc_curves.pdf"), width = 7, height = 4.9)
  corr.lists <- list()
  
  for (method_i in c(1:length(methods.notnulls))) {
    method.notnulls <- methods.notnulls[[method_i]]
    
    if (method.notnulls == "perseus") {
      old.method.name = "perseus"
    } else {
      old.method.name <- paste0(method.notnulls, ".", str_remove(method.allnulls, ".*\\.")) 
    }
    
    print(old.method.name)
    
    compare.to.use = (kOldImputationMethods[[old.method.name]]$sep.method == "mnar.1.nonzero")
    
    corr.lists[[method.notnulls]] <- 
      CountCorrRandomImputed(intensities = intensities,
                             sep.method = kOldImputationMethods[[old.method.name]]$sep.method,
                             method.mar = kOldImputationMethods[[old.method.name]]$method.mar,
                             method.mnar = kOldImputationMethods[[old.method.name]]$method.mnar,
                             use_only_not_nulls = T,
                             compare.with.mnar.1.nonzero = compare.to.use,
                             sigma_const = 4, metadata = metadata)
  }
  
  par(mar = c(5, 5, 3, 0))
  layout(matrix(1:4, nrow = 2), widths = c(.5, .5))
  
  colours <- unname(kColorList[c(4:(length(corr.lists) + 2), 2)])
  
  for (method_i in c(1:length(corr.lists))) {
    col = colours[[method_i]]
    MakeRocCurve(deltas.random = corr.lists[[method_i]]$corr.rand, 
                 deltas.imputed = corr.lists[[method_i]]$corr.imputed,
                 title = paste0(methods.notnulls[method_i], " (",
                                method.allnulls, ")"),
                 reverse = F, col = col, add = F, plot.hist = T)
  }
  
  aucs <- list()
  
  layout(matrix(1:2, nrow = 1), widths = c(.7, .3))
  
  for (method_i in c(1:length(corr.lists))) {
    if (method_i > 1) {
      add = T
    } else {
      add = F
    }
    
    col = colours[[method_i]]
    aucs[paste0(methods.notnulls[method_i])] <- 
      MakeRocCurve(deltas.random = 
                     corr.lists[[method_i]]$corr.rand, 
                   deltas.imputed = corr.lists[[method_i]]$corr.imputed,
                   title = "ROC for correlation",
                   xlab = "False positive rate", 
                   ylab = "True positive rate", 
                   reverse = F, col = col, add = add, plot.hist = F)
    aucs.df <- data.frame(auc = unlist(aucs))
  }
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = names(corr.lists),
         col = colours, 
         lty=1, bty = "n", lwd = 2)
  
  plot.new()
  aucs.df$auc <- format(aucs.df$auc, nsmall = 2)
  grid.table(aucs.df)
  
  dev.off()
}

CountAUCsAllMethods <- function(intensities, n.threads = 5, n.iterations = 2, 
                                use_only_not_nulls, sigma_const = 4, method.allnulls = 'method.1',
                                metadata) {
  
  aucs <- pbmclapply(c(1:n.iterations), function(seed.) {
    corr.lists <- list()
    
    for (method.notnulls in methods.notnulls) {
      
      if (method.notnulls == "perseus") {
        old.method.name = "perseus"
      } else {
        old.method.name <- paste0(method.notnulls, ".", str_remove(method.allnulls, ".*\\.")) 
      }
      
      print(old.method.name)
      
      compare.to.use = (kOldImputationMethods[[old.method.name]]$sep.method == "mnar.1.nonzero")
      
      corr.lists[[method.notnulls]] <- 
        CountCorrRandomImputed(intensities = intensities,
                               sep.method = kOldImputationMethods[[old.method.name]]$sep.method,
                               method.mar = kOldImputationMethods[[old.method.name]]$method.mar,
                               method.mnar = kOldImputationMethods[[old.method.name]]$method.mnar,
                               use_only_not_nulls = T,
                               compare.with.mnar.1.nonzero = compare.to.use,
                               sigma_const = 4, metadata = metadata)
    }
    
    corr.lists <- rapply(corr.lists, f = function(x) ifelse(is.na(x), 0 , x), how = "replace" )
    aucs.tmp <- list()
    
    for (method.not.nulls in methods.notnulls) {
      print(method.not.nulls)
      
      aucs.tmp[[method.not.nulls]] <- 
        MakeRocCurve(deltas.random = corr.lists[[method.not.nulls]]$corr.rand, 
                     deltas.imputed = corr.lists[[method.not.nulls]]$corr.imputed,
                     reverse = F, plot.hist = F)
    }    
    
    aucs.tmp
  }, mc.cores = n.threads, mc.silent = F)
  
  
  return(aucs)
}

MakeAUCsConfIntervals <- function(max.quant.obj, n.threads = 1, n.iterations = 2) {
  intensities <- max.quant.obj$GetIntensityTable()
  metadata <- max.quant.obj$GetMetadataExperiments()
  
  aucs <- CountAUCsAllMethods(intensities = intensities, n.threads = n.threads, 
                              n.iterations = n.iterations, 
                              use_only_not_nulls = T, sigma_const = 4, metadata = metadata)
  
  saveRDS(aucs, paste('./out/imputation/', MakeDate(), 'iBAQ_aucs.rds'))
  
  aucs.dt <- as.data.table(t(matrix(unlist(aucs), ncol = length(aucs), byrow = F)))
  colnames(aucs.dt) <- names(aucs[[1]])
  stat <- rbind('means' = colMeans(aucs.dt), 'sds' =  apply(aucs.dt, 2, sd)) 
  stat <- as.data.frame(t(stat))
  stat <- stat[order(stat$means, decreasing = T),]
  conf_intervals <- data.frame('method' = rownames(stat), 
                               'AUC' = paste(round(stat$means, 2), 
                               '±', round(qt(0.975, df = nrow(aucs.dt) - 1) *
                                            sqrt(stat$sds / nrow(aucs.dt)), 2), sep = ''))
  
  pdf(paste0('out/imputation/', MakeDate(), "_AUC_conf_intervals.pdf"), width = 5, height = 7)
  grid.table(conf_intervals, rows = NULL)
  dev.off()
  write.table(conf_intervals, paste('./out/imputation/', MakeDate(),
                                    'confidence_interval_AUC.tsv'), 
              row.names = F, sep = '\t', quote = F)
}



