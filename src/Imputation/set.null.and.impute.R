library(data.table)
library(stringr)
library(bit64)
library(scales)

source("lib/constants.R")

ImputeRandom <-
  function(dt,
           random.method,
           cur.exp.cols,
           vector.original.int) {
    
    tmp.dt = copy(dt)
    
    original.distribution <- vector.original.int[vector.original.int > 0]
    
    if (random.method == "rnorm") {
      mu.original <- mean(original.distribution)
      sd.original <- sd(original.distribution)
      for (col in cur.exp.cols) {
        data.table::set(
          tmp.dt,
          i = which(tmp.dt[[col]] == 0),
          j = col,
          value = rnorm(length(which(tmp.dt[[col]] == 0)),
                        mu.original, sd.original)
        )
      }
    } else if (random.method == "runif") {
      for (col in cur.exp.cols) {
        data.table::set(
          tmp.dt,
          i = which(tmp.dt[[col]] == 0),
          j = col,
          value = runif(
            length(which(tmp.dt[[col]] == 0)),
            min = min(original.distribution, na.rm = T),
            max = max(original.distribution, na.rm = T)
          )
        )
      }
    }
    
    tmp.dt
  }

MakeRandomAndImputedMtx <- function(intensities,
                                    cur.exp.cols, sep.method,
                                    method.mar, method.mnar,
                                    sigma.const = 4,
                                    random.method = "rnorm",
                                    part.to.delete = 0.5, seed,
                                    compare.with.mnar.1.nonzero = F) {
  
  if (!missing(seed)) 
    set.seed(seed) 
                                      
  cur.dt.original <- intensities[, c(cur.exp.cols, "Protein IDs"), with = F]

  cur.dt.original <- cur.dt.original[(rowSums(cur.dt.original[, cur.exp.cols, with = F] == 0) <
                                        dim(cur.dt.original)[[2]] - 1),]
  
  vector.original.int <- unlist(cur.dt.original[, cur.exp.cols, with = F])
  cur.dt.original.deleted = copy(cur.dt.original)
  

  rows.deleted <- c()
  cols.deleted <- c()
  
  for (col.count in c(1:dim(cur.dt.original[, cur.exp.cols, with = F])[[2]])) {
    non.zero.ind <- which(cur.dt.original[, col.count, with = F] > 0)
    count.to.delete <- floor(length(non.zero.ind) * part.to.delete)
    rows.deleted <-
      c(rows.deleted,
        sample(non.zero.ind, count.to.delete, replace = F))
    cols.deleted <-
      c(cols.deleted, rep(col.count, count.to.delete))
  }
  
  
  for (i in c(1:length(rows.deleted))) {
    cur.dt.original.deleted[rows.deleted[i], cols.deleted[i]] <- 0
  }
  
  if (sep.method == 'mnar.1.nonzero' | compare.with.mnar.1.nonzero) {
    rows.all.nulls <-
      c(cur.dt.original.deleted[rowSums(cur.dt.original.deleted == 0) == 
                                  (ncol(cur.dt.original.deleted) - 1), which = T], 
        cur.dt.original.deleted[rowSums(cur.dt.original.deleted == 0) == 
                                  (ncol(cur.dt.original.deleted) - 2), which = T])
  } else {
    rows.all.nulls <-
      cur.dt.original.deleted[rowSums(cur.dt.original.deleted == 0) == 
                                (ncol(cur.dt.original.deleted) - 1), which = T] 
  }
    
  cur.dt.random <-
    ImputeRandom(cur.dt.original.deleted,
                 "runif",
                 cur.exp.cols,
                 vector.original.int)
  
  cur.dt.imputed <-
    Impute(cur.dt.original.deleted,
           'Protein IDs',
           sep.method = sep.method,
           method.mar = method.mar,
           method.mnar = method.mnar,
           sigma.const = sigma.const)
  
  # cur.dt.imputed to normal row order
  
  old.row.order <- match(cur.dt.original$`Protein IDs`, cur.dt.imputed$`Protein IDs`)
  cur.dt.imputed <- cur.dt.imputed[old.row.order, ]
  identical(cur.dt.imputed$`Protein IDs`, cur.dt.original$`Protein IDs`)
  
  # make two distributions of deltas
  
  mtx.original <- as.matrix(cur.dt.original[, cur.exp.cols, with = F])
  mtx.random <- as.matrix(cur.dt.random[, cur.exp.cols, with = F])
  mtx.imputed <- as.matrix(cur.dt.imputed[, cur.exp.cols, with = F])
  
  mtx.list <-
    list(
      original = mtx.original,
      random = mtx.random,
      imputed = mtx.imputed,
      rows.deleted = rows.deleted,
      cols.deleted = cols.deleted,
      rows.all.nulls = rows.all.nulls
    )
  
  mtx.list
}


