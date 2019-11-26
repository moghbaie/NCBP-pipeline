require(data.table)
require(digest)
require(dendextend)
require(stringr)
stringsAsFactors = FALSE

SetNodesColor <- function(node, metadata) {
  if (!is.leaf(node)) {
   return(node) 
  }
  
  if (is.na(attributes(node)$label)) {
    color <- kColorList[["no group"]]
  } else {
    color <- 
      kColorList[[unlist(metadata[metadata$Protein.IDs == attributes(node)$label, "group"])]]
  }
  
  attr(node, "nodePar") <- 
      c(attributes(node)$nodePar,list(cex=1.5, lab.cex=0.5, pch=20, col=color, lab.font=1))
  
  return(node)
}

MakeDendroWithClusters = function(dist, metadata, n.of.clusters = 5) {
  int.hc <- hclust(dist)
  # save cluster for behaviour plot
  clusters <- cutree(int.hc, k = n.of.clusters)
  dend <- as.dendrogram(int.hc)
  # save order of clusters for behaviour plot
  labels.order <- rev(labels(dend))
  clusters.order <- unique(unname(clusters[labels.order]))
  # set colors to nodes according to groups
  dend <- dendrapply(dend, function(node) {
    SetNodesColor(node, metadata)
  })
  # convert protein ids to gene names
  gene.names <- metadata[match(labels(dend), metadata$Protein.IDs), Gene.names]
  dend <- dend %>% set("labels", gene.names)  
  result <- list('dend' = dend, 'clusters_order' = clusters.order, 'clusters' = clusters)
  return(result)
} 
