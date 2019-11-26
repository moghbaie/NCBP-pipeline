source("src/Clusterization/make_dendrogram.R")

MaxQuant$set("public", "MakeClusters", function(clustering.method = "hclust", n.of.clusters = 9) {
  if (clustering.method == "hclust") {
    private$clusters <- MakeDendroWithClusters(private$distance.matrix, private$metadata.proteins,
                                               n.of.clusters)$clusters
  } else {
    stop("Unknown clustering method")
  }
  
  private$AddChange(paste0("clustered with method ", clustering.method, 
                           ". Number of clusters: ", n.of.clusters))
  invisible(self)
})
