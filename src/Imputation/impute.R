library(data.table)
library(stringr)
library(bit64)
library(reshape2)

source("./lib/constants.R")

MergeAll <- function(x, y)
  merge(x, y, all = TRUE)


CountMeanSDNoOutliers <- function(sample) {
  OutVals = boxplot(sample, plot = F)$out
  sample <- sample[!sample %in% OutVals]
  mu <- mean(sample, na.rm=T)
  sd <- sd(sample, na.rm=T)
  
  list("mu" = mu, "sd" = sd)
}


ImputeMAR <- function(df.input, id.column, method, sigma.const, stat.info) {
  
  print("Impute MAR")
  
  if (empty(df.input)) {
    return(df.input)
  }
  
  df <- copy(df.input[, -id.column, with=F])
  
  df[df == 0] <- NA
  
  meta.imputation <- data.table(id = rownames(df), 
                                mean.intensity = rowMeans(df, na.rm=T),
                                n.not.zero = rowSums(!is.na(df)),
                                n.zero = rowSums(is.na(df)))
  
  if (method %in% c("mar.stat.1", "mar.stat.2", "mar.stat.3", "mar.stat.4", "mar.stat.5", 
                    "mar.stat.6")) {
    for (row in c(1:nrow(df))) {
      for (column in colnames(df)) {
        
        if (is.na(df[row, column, with = F])) {
          corrs <- c()
          
          for (other.column in colnames(df)) {
            if (!is.na(df[row, other.column, with = F])) {
              corrs <- c(corrs, stat.info$cor.mtx[column, other.column])
            }
          }
          
          corrs.mean <-
            if (abs(mean(corrs)) == 0)
              0.1
            else
              mean(corrs)
          
           if (method == "mar.stat.1") {
            delta.new <-
                rnorm(1, mean = stat.info$mu.deltas, sd = sqrt(2) * stat.info$sd.deltas /
                        abs(corrs.mean))
            
            int.new <-
              meta.imputation[row]$mean.intensity * abs(1 + delta.new)
            
          } else if (method == "mar.stat.2") {
            delta.new <-
              rnorm(
                1,
                mean = stat.info$mu.deltas,
                sd = sqrt(2) * stat.info$sd.deltas / abs(corrs.mean) * meta.imputation[row]$n.zero
              )
            
            int.new <-
              meta.imputation[row]$mean.intensity * abs(1 + delta.new)
            
          } else if (method == "mar.stat.3") {
            delta.new <-
              rnorm(1,
                    mean = stat.info$mu.deltas,
                    sd = sqrt(2) * stat.info$sd.deltas / abs(corrs.mean))
            
            int.new <-
              meta.imputation[row]$mean.intensity * abs(1 + delta.new) * 
              (1 - log10(meta.imputation[row]$n.zero))
            
          } else if (method == "mar.stat.4") {
            delta.new <-
              rnorm(1,
                    mean = stat.info$mu.deltas,
                    sd = sqrt(2) * stat.info$sd.deltas / abs(corrs.mean))
            
            int.new <-
              meta.imputation[row]$mean.intensity * abs(1 + delta.new) * 3 * 
              (meta.imputation[row]$n.not.zero) /
              ncol(df)
            
          } else if (method == "mar.stat.5") {
            delta.new <-
              rnorm(
                1,
                mean = stat.info$mu.deltas,
                sd = sqrt(2) * stat.info$sd.deltas / abs(corrs.mean) * 
                  (1 + log10(meta.imputation[row]$n.zero))
              )
            
            int.new <-
              meta.imputation[row]$mean.intensity * abs(1 + delta.new)
            
          } else if (method == "mar.stat.6") {
            delta.new <-
              rnorm(1,
                    mean = stat.info$mu.deltas,
                    sd = (1 / sqrt(2)) * stat.info$sd.deltas / abs(corrs.mean))
            
            int.new <-
              meta.imputation[row]$mean.intensity * abs(1 + delta.new)
          }
          
          data.table::set(df, row, column, int.new)
        }
      }
    }
    
  } else {
    stop(paste("Unknown method:", method))
  }
  
  df[[id.column]] = df.input[[id.column]]
  
  return(df)
}


ImputeMNAR <- function(df.input, id.column, method, sigma.const, stat.info, colnames.to.take) {
  
  print("Impute MNAR")
  
  if (empty(df.input)) {
    return(df.input)
  }
  
  df <- copy(df.input[, -id.column, with=F])
  
  if (missing(colnames.to.take)) {
    colnames.to.take <- colnames(df)
  }

  # different imputation methods
  
  if (method == "perseus" | method == "mnar.stat.2") {
    for (colname in colnames.to.take) {
      column <- unlist(df[, colname, with = F])
      
      colname.bool <- rownames(stat.info$column.stats) == colname
      
      data.table::set(df, which(column == 0), colname,
                      rnorm(
                        length(which(column == 0)),
                        mean = stat.info$column.stats[colname.bool, ]$mu - 
                          1.8 * stat.info$column.stats[colname.bool, ]$sd,
                        sd = 0.3 * stat.info$column.stats[colname.bool, ]$sd
                      ))
    }
  } else if (method == "mnar.stat.1") {
    
    for (colname in colnames.to.take) {
      column <- unlist(df[, colname, with = F])
      
      colname.bool <- rownames(stat.info$column.stats) == colname
      
      data.table::set(df, which(column == 0), colname,
                      runif(
                        length(which(column == 0)),
                        1, 
                        stat.info$column.stats[colname.bool, ]$mu - 
                          sigma.const * stat.info$column.stats[colname.bool, ]$sd
                      ))
    }
  } else {
    stop(paste("Unknown method:", method))
  }
  
  df[[id.column]] = df.input[[id.column]]
  
  return(df)
}


SeparateByNNonZero <- function(data, id.column, n.nonzero = 0) {
  is.mnar <- rowSums(data[, -id.column, with=F] > 0) <= n.nonzero
  
  list(mar = data[!is.mnar, ], mnar = data[is.mnar, ])
}


SeparateImputeByNNonZeroMerge <- function(data, id.column, n.nonzero = 0, stat.info,
                                          method.mar = "mar.stat.6",
                                          method.mnar = 'mnar.stat.1',
                                          sigma.const = 4, remove.outliers = F) {
  
  full.rows <- rowSums(data[, -id.column, with = F] == 0) == 0
  
  sep.data <- SeparateByNNonZero(data[!full.rows, ], id.column, n.nonzero)
  
  imputed.mar <- ImputeMAR(sep.data$mar, id.column, method.mar, sigma.const, stat.info)
  imputed.mnar <- ImputeMNAR(sep.data$mnar, id.column, method.mnar, sigma.const, stat.info)
  
  return(rbindlist(list(imputed.mar, imputed.mnar, data[full.rows, ])))
}


SeparateImputeOneColForMNAR <- function(df, id.column, n.nonzero, stat.info, 
                                        method.mar, method.mnar,
                                        sigma.const, remove.outliers) {
  
  n.not.null.proteins <- colSums(df[, -id.column, with = F] != 0)
  
  full.rows <- rowSums(df[, -id.column, with = F] == 0) == 0
  sep.data <- SeparateByNNonZero(df[!full.rows, ], id.column, n.nonzero)
  
  imputed.mnar <- ImputeMNAR(sep.data$mnar, id.column, method.mnar, sigma.const, stat.info,
                             colnames.to.take = names(which.max(n.not.null.proteins)))
  imputed <- ImputeMAR(rbindlist(list(sep.data$mar, imputed.mnar)), 
                       id.column, method.mar, sigma.const, stat.info)
  
  return(rbindlist(list(imputed, df[full.rows, ])))
}


SeparateImputeMARMNARMerge <- function(df, id.column, sep.method = "mnar.1.nonzero",
                                       method.mar = "mar.stat.6",
                                       method.mnar = 'mnar.stat.1',
                                       sigma.const = 4, remove.outliers = F) {

  # calculate correlation using only pairwise elements
  
  df.with.na <- copy(df[, -id.column, with = F])
  df.with.na[df.with.na == 0] <- NA 
  
  cor.mtx = cor(df.with.na, use = "pairwise.complete.obs", method = "pearson")
  cor.mtx[is.na(cor.mtx)] <- 0
  
  # calculate delta distribution
  
  epairs <- t(combn(colnames(df[, -id.column, with = F]), m = 2))
  deltas <- c()
  
  for (pair.id in 1:nrow(epairs)) {
    pair.deltas <-
      (df.with.na[, epairs[pair.id, 1], with = F] - df.with.na[, epairs[pair.id, 2], with = F]) /
      rowMeans(df.with.na[, epairs[pair.id,], with = F])
    pair.deltas <- unlist(pair.deltas)
    deltas <- c(deltas, pair.deltas[!is.na(pair.deltas)])
  }
  
  mu.delta <- mean(deltas)
  sd.delta <- sd(deltas)
  
  if (is.na(sd.delta)) {
    sd.delta <- 0
  }
  
  if (remove.outliers) {
    column.stats <-
      rbindlist(apply(df[, -id.column, with = F], 2, function(column) {
        CountMeanSDNoOutliers(column[column != 0])
      }))
  } else {
    column.stats <-
      rbindlist(apply(df[, -id.column, with = F], 2, function(column) {
        list("mu" = mean(column[column != 0], na.rm=T), "sd" = sd(column[column != 0], na.rm = T))
      }))
  }
  rownames(column.stats) <- colnames(df[, -id.column, with = F])
  
  # make stat.info
  stat.info <- list(mu.deltas = mu.delta, sd.deltas = sd.delta, cor.mtx = cor.mtx, 
                    column.stats = column.stats)
  
  # different separation methods
  if (sep.method == "all.mnar") {
    return(ImputeMNAR(df, id.column, method = method.mnar, sigma.const = sigma.const, 
                      stat.info = stat.info))
  } else if (sep.method == "mnar.0.nonzero") {
    return(SeparateImputeByNNonZeroMerge(df, id.column, 0, stat.info = stat.info, 
                                         method.mar = method.mar,
                                         method.mnar = method.mnar, sigma.const = sigma.const,
                                         remove.outliers = remove.outliers))
  } else if (sep.method == "mnar.1.nonzero") {
    return(SeparateImputeByNNonZeroMerge(df, id.column, 1, stat.info = stat.info, 
                                         method.mar = method.mar,
                                         method.mnar = method.mnar, sigma.const = sigma.const,
                                         remove.outliers = remove.outliers))
  } else if (sep.method == "impute.1.col.for.mnar.0.nonzero") {
    return(SeparateImputeOneColForMNAR(df, id.column, 0, stat.info = stat.info, 
                                       method.mar = method.mar,
                                       method.mnar = method.mnar, sigma.const = sigma.const,
                                       remove.outliers = remove.outliers))
  } else if (sep.method == "impute.1.col.for.mnar.1.nonzero") {
    return(SeparateImputeOneColForMNAR(df, id.column, 1, stat.info = stat.info, 
                                       method.mar = method.mar,
                                       method.mnar = method.mnar, sigma.const = sigma.const,
                                       remove.outliers = remove.outliers))
  } else {
    stop(paste("Unknown sep.method:", sep.method))
  }

}

# impute 1 case/control from 1 experiment
Impute <- function(df.input, id.column = NA,
                   sep.method = "mnar.1.nonzero",
                   method.mar = 'mar.stat.6',
                   method.mnar = "mnar.stat.1",
                   sigma.const = 4, remove.outliers = F) {
    
    df <- copy(df.input)
    
    df.imputed <-
      SeparateImputeMARMNARMerge(df = df, id.column = id.column, 
                                 sep.method = sep.method, method.mar = method.mar,
                                 method.mnar = method.mnar, sigma.const = sigma.const,
                                 remove.outliers = remove.outliers)
    
    return(df.imputed)
}


# impute case AND control from 1 group
ImputeCaseAndControl <-
  function(intensities.group,
           experiments.case,
           experiments.control,
           id.column = NA,
           sep.method = "mnar.1.nonzero",
           method.mar = "mar.stat.6",
           method.mnar = 'mnar.stat.1',
           sigma.const = 4, remove.outliers = F) {
    
    intensities.group <-
      intensities.group[rowSums(intensities.group[, c(experiments.case, experiments.control), 
                                                  with = F]) != 0, ]
    
    intensities.group.case <-
        intensities.group[, c(experiments.case, id.column), with = F]
      
    intensities.group.control <-
        intensities.group[, c(experiments.control, id.column), with = F]
    
    imputed.data <- lapply(list(intensities.group.case, intensities.group.control),
                           function(df) {
                             Impute(
                                 df.input = df,
                                 id.column = id.column,
                                 sep.method = sep.method,
                                 method.mar = method.mar,
                                 method.mnar = method.mnar,
                                 sigma.const = sigma.const
                               )
                           })
    
    print('merge')
    
    intensities.group.imputed <- merge(
      imputed.data[[1]],
      imputed.data[[2]],
      by = id.column,
      all = T
    )
    
    return(intensities.group.imputed)
  }


# Imputation for the whole data set
MaxQuant$set("public", "ImputeAll", function(id.column = NA,
                                             sep.method = "mnar.1.nonzero",
                                             method.mar = "mar.stat.6",
                                             method.mnar = 'mnar.stat.1',
                                             sigma.const = 4,
                                             save.result.table = T,
                                             path. = 'out/', 
                                             file.name, seed, remove.outliers=F) {
  if (!missing(seed)) 
    set.seed(seed) 
  
  if (is.na(id.column)) {
    stop('ImputeAll needs id.column to perform imputation')
  }
  
  private$is.zero.intensity <- private$intensity.table[, -kIdColumnName, with = F]
  private$is.zero.intensity[, names(private$is.zero.intensity) := lapply(.SD, as.logical)]
  private$is.zero.intensity[[kIdColumnName]] <- private$intensity.table[[kIdColumnName]]
  
  
  # TODO: move it to BuildMaxQuant
  private$metadata.experiments$group <-
    paste(
      private$metadata.experiments$location,
      private$metadata.experiments$target,
      private$metadata.experiments$buffer,
      sep = '_'
    )
  
  out.list <-
    lapply(unique(private$metadata.experiments$group), function(exp.group) {
      print(exp.group)
      columns.to.take <-
        c(private$metadata.experiments[private$metadata.experiments$group == exp.group, 
                               ]$experiment.name, id.column)
      
      intensities.group <-
        private$intensity.table[, columns.to.take, with = F]
      
      ImputeCaseAndControl(
        intensities.group = intensities.group,
        experiments.case = 
          private$metadata.experiments[private$metadata.experiments$group == exp.group &
                                      private$metadata.experiments$is.case == T, ]$experiment.name,
        experiments.control = 
          private$metadata.experiments[private$metadata.experiments$group == exp.group &
                                       private$metadata.experiments$is.case == F, ]$experiment.name,
        id.column = id.column,
        sep.method = sep.method,
        method.mar = method.mar,
        method.mnar = method.mnar,
        sigma.const = sigma.const, remove.outliers = remove.outliers
      )
    })
  
  private$intensity.table <- Reduce(MergeAll, out.list)
  
  private$AddChange("imputed")
  
  if (save.result.table) {
    date <- MakeDate()
    
    if (missing(file.name)) {
      file.name = paste(sep.method, method.mar, method.mnar, sep = "_") 
    }
    
    saveRDS(
      private$intensity.table,
      paste0(
        path.,
        date,
        '_imputed_',
        file.name,
        '.rds'
      )
    )
  }
  invisible(self)
})
