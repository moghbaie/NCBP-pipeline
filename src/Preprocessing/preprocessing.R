require(stringr)
require(scales)
require(data.table)
source('lib/constants.R')
source('lib/helper_functions.R')

MaxQuant$set("public", "SumMAGOH", function(save.result.table = T) {
  sum.values <- 
    private$intensity.table[grep('P61326', private$intensity.table$`Protein IDs`),
                            private$metadata.experiments$experiment.name, with = F] + 
    private$intensity.table[grep('Q96A72', private$intensity.table$`Protein IDs`),
                            private$metadata.experiments$experiment.name, with = F] 
  
  data.table::set(private$intensity.table, 
                  grep('P61326', private$intensity.table$`Protein IDs`),
                  private$metadata.experiments$experiment.name, as.list(sum.values)) 
  
  private$metadata.proteins[grep('P61326', 
                                 private$metadata.proteins$Protein.IDs),]$Razor...unique.peptides <- 
    private$metadata.proteins[grep('P61326', 
                                   private$metadata.proteins$Protein.IDs),]$Razor...unique.peptides + 
    private$metadata.proteins[grep('Q96A72', 
                                   private$metadata.proteins$Protein.IDs),]$Razor...unique.peptides
  
  private$intensity.table <- 
    private$intensity.table[!grepl('Q96A72', private$intensity.table$`Protein IDs`),]
  private$metadata.proteins <- 
    private$metadata.proteins[!grepl('Q96A72', private$metadata.proteins$Protein.IDs),]
  
  private$AddChange("summed MAGOH and MAGOHB")
  
  if (save.result.table) {
    date <- paste(str_split_fixed(Sys.time(), ' ', 3)[,1], str_split_fixed(Sys.time(), ' ', 3)[,2],
                  sep='_')
    saveRDS(private$intensity.table, paste0('out/', date, '_magoh_sumed',  '.rds')) 
  }
  invisible(self)
})

MaxQuant$set("public", "RemoveContaminantsReversed", function(intensities, metadata) {
  good.ids <- 
    private$metadata.proteins[(private$metadata.proteins$contaminant.bool == F)& 
                                (private$metadata.proteins$reverse.bool == F),]$Protein.IDs
  private$intensity.table <-
    private$intensity.table[private$intensity.table$`Protein IDs` %in% good.ids,]
  
  private$metadata.proteins <-
    private$metadata.proteins[private$metadata.proteins$Protein.IDs %in% good.ids, ]
  
  private$AddChange("removed potential comtaminants and reversed")
  invisible(self)
})

MaxQuant$set("public", "LogTricky", function() {
  numeric.columns <- unlist(lapply(private$intensity.table, is.numeric))
  numeric.columns <- names(numeric.columns[numeric.columns == TRUE])
  
  intensities.tmp <- log2(private$intensity.table[, numeric.columns, with = F])
  
  for (j in 1:ncol(intensities.tmp)) {
    data.table::set(intensities.tmp, which(is.infinite(intensities.tmp[[j]])), j, 0)
  }
  
  protein.id.column <- private$intensity.table[[kIdColumnName]]
  private$intensity.table <- intensities.tmp
  private$intensity.table[[kIdColumnName]] <- protein.id.column

  private$AddChange("log2")
  invisible(self)
})

MaxQuant$set("public", "PowTwo", function() {
  numeric.columns <- unlist(lapply(private$intensity.table, is.numeric))
  numeric.columns <- names(numeric.columns[numeric.columns == TRUE])
  
  intensities.tmp <- as.data.table(2^(private$intensity.table[, numeric.columns, with = F]))
  
  protein.id.column <- private$intensity.table[[kIdColumnName]]
  private$intensity.table <- intensities.tmp
  private$intensity.table[[kIdColumnName]] <- protein.id.column
  
  private$AddChange("2^{intensity}")
  invisible(self)
})

MaxQuant$set("public", "Normalize", function(protein.id = 'Q9U6Y5',
                                             norm.method = 'by_protein', 
                                             is.log.transformed = F) {
  
    numeric.columns <- unlist(lapply(private$intensity.table, is.numeric))
    numeric.columns <- names(numeric.columns[numeric.columns == TRUE])
    
    intensity.tmp <- data.table()
    
    if (is.log.transformed) {
      if (norm_method == 'by_protein') {
        intensity.tmp <- 
          apply(private$intensity.table[, numeric.columns, with = F], 2, function(col) {
              col - col[grep(protein.id, unlist(private$intensity.table[, kIdColumnName, with = F]))]
        })
      } else {
        stop("Unknown normalization method")
      }
    } else {
      if (norm.method == 'by_protein') {
        intensity.tmp <- 
          apply(private$intensity.table[, numeric.columns, with = F], 2, function(col) {
            col/col[grep(protein.id, unlist(private$intensity.table[, kIdColumnName, with = F]))]
        })
      } else if(norm.method == 'scale'){
        
        intensity.tmp <- private$intensity.table[, numeric.columns, with = F] / 1
        for (j in 1:ncol(intensity.tmp)) {
          data.table::set(intensity.tmp, which(intensity.tmp[[j]] == 0), j, value = NA)
        }  
        intensity.tmp <- as.data.table(scale(as.matrix(log10(intensity.tmp)))) + 4
        intensity.tmp <- 10 ^ intensity.tmp
        for (j in 1:ncol(intensity.tmp)) {
          data.table::set(intensity.tmp, which(is.na(intensity.tmp[[j]])), j, value = 0)
        }  
      }
      else {
        stop("Unknown normalization method")
      }
    }  
    
    intensity.tmp <- cbind(intensity.tmp, private$intensity.table[, kIdColumnName, with = F])
    private$intensity.table <- intensity.tmp
    
    change.str <- paste0("normalized with method ", norm.method)
    if (norm.method == 'by_protein') {
      change.str <- paste0(change.str, " (protein id: ", protein.id, ")")
    }
    private$AddChange(change.str)
    
    invisible(self)
})

MaxQuant$set("public", "Rescale", function(range = c(10, 1000), by.coef = F) {
  
  if (length(range) != 2) {
    stop("range length must be 2")
  }
  
  numeric.columns <- unlist(lapply(private$intensity.table, is.numeric))
  numeric.columns <- names(numeric.columns[numeric.columns == TRUE])
  intensities.m <- as.matrix(private$intensity.table[, numeric.columns, with = F])
  
  if (by.coef) {
    m <- mean(rowMeans(private$intensity.table[, numeric.columns, with = F], na.rm = T))
    k <- mean(range) / m
    intensities.m <- intensities.m * k
    private$AddChange(paste('rescaled by coefficient ', k, sep = ''))
  } else{
    intensities.m <- rescale(x = intensities.m, to = range)  
    private$AddChange(paste("rescaled to range", paste(range, collapse = " ")))
  }
  protein.id.column <- private$intensity.table[[kIdColumnName]]
  private$intensity.table <- as.data.table(intensities.m)
  private$intensity.table[[kIdColumnName]] <- protein.id.column
  
  invisible(self)
})

FindIntensityRange <- function(max.quant.obj) {
  numeric.columns <- unlist(lapply(max.quant.obj$GetIntensityTable(), is.numeric))
  numeric.columns <- names(numeric.columns[numeric.columns == TRUE])
  
  return(c(min(max.quant.obj$GetIntensityTable()[, numeric.columns, with = F], na.rm = T),
           max(max.quant.obj$GetIntensityTable()[, numeric.columns, with = F], na.rm = T)))
}

MaxQuant$set("public", "FilterColumns", function(pattern) {
  cols.to.drop = colnames(private$intensity.table)[grep(pattern, colnames(private$intensity.table))]
  private$intensity.table <- private$intensity.table[, -cols.to.drop, with = F]
  
  if (!is.null(private$pvals.logfc)) {
    cols.to.drop = colnames(private$pvals.logfc)[grep(pattern, colnames(private$pvals.logfc))]
    private$pvals.logfc <- private$pvals.logfc[, -cols.to.drop, with = F]
  }

  rows.in.metadata <- grep(pattern, private$metadata.experiments$experiment.name)
  private$metadata.experiments <- private$metadata.experiments[-rows.in.metadata, ]
  
  private$AddChange(paste("deleted experiments with pattern", pattern))
  
  invisible(self)
})

CountMeanBetweenReplicasForCases <- function(intensity.dt, id.col = "Protein IDs", metadata,
                                             keep.id.column = F) {
  df <- data.frame(Protein.IDs = intensity.dt[, id.col, with=F])
  
  for (group. in unique(metadata$group)) {
    df[, paste0(group.)] <-
      rowMeans(subset(intensity.dt, select = 
                        metadata[metadata$group == group. & metadata$is.case == T, ]$experiment.name), 
               na.rm = TRUE)
  }
  
  rownames(df) <- df$Protein.IDs
  if (keep.id.column) {
    return(df)
  }
  
  return(df[, !(colnames(df) == "Protein.IDs")])
}

MakeDataFrameWithCases <- function(intensity.dt, metadata.experiments, id.col = "",
                                   keep.id.column = F) {
  cols <- unlist(metadata.experiments[(is.case), experiment.name])
    
  df <- as.data.frame(intensity.dt[, cols, with = F])
  df$Protein.IDs = as.vector(unlist(intensity.dt[, id.col, with=F]))
  rownames(df) <- df$Protein.IDs
  if (keep.id.column) {
    return(df)
  }
  
  return(df[, !(colnames(df) == "Protein.IDs")])
}

MaxQuant$set("public", "MakePrettyExperimentNames", function() {
  pretty.names <- fread('data/pretty_experiments_names.csv')
  metadata <- private$metadata.experiments
  metadata$experiment <- paste(metadata$location, metadata$target, sep = '_')
  metadata <- merge(metadata, pretty.names, by = 'experiment')
  metadata <- metadata[, names(metadata) %nin% c('experiment'), with = F]
  private$metadata.experiments <- metadata
  invisible(self)
})

