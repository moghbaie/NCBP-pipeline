require(proxy)
require(data.table)

MakeEuclideanDist <- function(df.num) {
  mtx <- as.matrix(df.num[, apply(df.num, 2, function(x) var(x, na.rm = TRUE)) != 0])
  return(dist(mtx))
}

MakeDistELife <- function(dt.num) {
  # mtx <- as.matrix(dt.num[, apply(dt.num, 2, var, na.rm = TRUE) != 0])
  mtx <- as.matrix(dt.num)
  d1.all <- dist(mtx, method = 'cosine')
  d1.all[d1.all < 0] <- 0
  d1.all[is.na(d1.all)] <- 1   # is it correct?
  d2.all <- dist(mtx, method = 'euclidean')
  d2.all <-
    scales::rescale(d2.all, to = c(0, 0.9))   # why magic 0.9?
  d2.all[is.na(d2.all)] <- 1   # is it correct?
  d12.all <- log(d1.all * d2.all)   # why log?
  d12.all <- scales::rescale(d12.all, to = c(0, 1))
  d12.all[is.infinite(d12.all)] <- 0
  return(as.dist(d12.all))
}

MakeDistAllTargets <- function(dt, id.col = "", targets = c(), metrics = "euclidean",
                                  all.targets = T) {
  
  df <- as.data.frame(dt) 
  if (id.col != "") {
    rownames(df) = as.vector(unlist(df[, id.col]))
  }
  
  if (all.targets) {
    targets = c(kAllTargetsPattern)
  }
  
  result <- list()
  
  for (target in targets) {
    df.target <- df[, grepl(target, colnames(df))]
    df.target <- df.target[rowSums(!is.na(df.target)) > 0, ]
    df.target[is.na(df.target)] <- 0
    if (metrics == "euclidean") {
      result[[target]] <- MakeEuclideanDist(df.target)
    } else if (metrics == "eLife") {
      result[[target]] <- MakeDistELife(df.target)
    } else {
      stop("Unknown metrics")
    }
  }
  
  return(result)
}

MaxQuant$set("public", "CountDistanceMatrix", function(metrics = "euclidean", all.targets = T,
                                                       targets = c()) {
  if (all.targets) {
    private$distance.matrix <- MakeDistAllTargets(private$intensity.table, id.col = kIdColumnName,
                                                  metrics = metrics, all.targets = all.targets)[[1]]
  } else {
    stop("Not implemented")
  }
  
  private$AddChange(paste0("Calculated distance matrix with metrics ", metrics, 
                            ", all.targets = ", all.targets))
  invisible(self)
})


