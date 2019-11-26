require(plotly)
require(ape)
require(vegan)
require(ggplot2)
source('lib/ggplot.colours.R')

PlotPCOA <- function(pcoa.df, x, y, title = '', subtitle = '', color_list = kColorList) {

  g <- ggplot(pcoa.df, aes_string(x = x, y = y)) +
    geom_point(
      data = filter(pcoa.df, group == "no group"),
      alpha = 0.5,
      color = "grey"
    ) +
    geom_point(data = filter(pcoa.df, group != "no group"),
               alpha = 1,
               aes(color = group)) +
    labs(title = title, subtitle = subtitle) + theme(aspect.ratio = 1) +
    scale_colour_manual(values = color_list) +
    geom_text_repel(
      data = filter(pcoa.df, group != "no group"),
      aes(label = gene.names),
      hjust = 0.5,
      vjust = 0,
      alpha = 0.5,
      segment.alpha = 0.2,
      size = 2
    ) + theme_classic()

  g
}

PlotMDSWithArrows <- function(max.quant.obj) {
  df.intensities <- MakeDataFrameWithCases(max.quant.obj$GetIntensityTable(),
                                           max.quant.obj$GetMetadataExperiments(),
                                           id.col = kIdColumnName, keep.id.column = F)
  dist.mtx <- max.quant.obj$GetDistanceMatrix()
  metadata <- max.quant.obj$GetMetadataProteins()
  metadata[is.na(metadata$uniprot.gene.best)]$uniprot.gene.best <- metadata[is.na(metadata$uniprot.gene.best)]$Gene.names
  metadata[is.na(metadata$uniprot.gene.best)]$uniprot.gene.best <- metadata[is.na(metadata$uniprot.gene.best)]$Protein.IDs
  metadata$uniprot.gene.best <- str_split_fixed(metadata$uniprot.gene.best, ';', 2)[, 1]

  pcoa <- cmdscale(dist.mtx, k = 3, add = F)

  efit <- envfit(pcoa, df.intensities, choices = c(1:3))
  pcoa.df <- as.data.frame(pcoa)
  colnames(pcoa.df) <- c('Axis.1','Axis.2', 'Axis.3')
  rows.in.metadata <-
    match(rownames(pcoa.df), metadata$Protein.IDs)
  rows.in.metadata <- rows.in.metadata[!is.na(rows.in.metadata)]
  pcoa.df$group = factor(as.character(metadata[rows.in.metadata, ][["group"]]))
  pcoa.df$gene.names = as.character(metadata[rows.in.metadata, ][["uniprot.gene.best"]])
  pcoa.df$Protein.IDs <- rownames(pcoa.df)

  arrows <- as.data.frame(efit$vectors$arrows)
  min.arrow = max(min(pcoa.df$Axis.1), min(pcoa.df$Axis.2), min(pcoa.df$Axis.3))
  max.arrow = min(max(pcoa.df$Axis.1), max(pcoa.df$Axis.2), max(pcoa.df$Axis.3))
  arrow.range <- c(min.arrow, max.arrow)
  
  arrows$Dim1 <- rescale(arrows$Dim1, arrow.range) / 2
  arrows$Dim2 <- rescale(arrows$Dim2, arrow.range) / 2
  arrows$Dim3 <- rescale(arrows$Dim3, arrow.range) / 2
  arrows$experiments <- rownames(arrows)
  colnames(arrows) <- c('Axis.1', 'Axis.2', 'Axis.3', 'experiments')
  arrows$target <- as.factor(str_split_fixed(arrows$experiments,'_',3)[,2])
  arrows <- arrows[grep('_C_', arrows$experiments, invert = T), ]
  
  metadata.experiments <- max.quant.obj$GetMetadataExperiments()  
  arrows$experiments <- paste(metadata.experiments[match(arrows$experiments, 
                    metadata.experiments$experiment.name), ]$pretty_name, str_split_fixed(arrows$experiments, '_', 2)[, 2], sep = '_')
  
  pdf(paste0('out/', MakeDate(), '_mds_arrows.pdf'), width = 10, height = 10)

  axis.combinations <- combn(c('Axis.1', 'Axis.2', 'Axis.3'), 2)
  for (col.count in c(1:ncol(axis.combinations))) {
    current.x <- axis.combinations[1, col.count]
    current.y <- axis.combinations[2, col.count]

    g <- PlotPCOA(pcoa.df, current.x, current.y) +
      new_scale_color() +
      geom_segment(data = arrows,
                   aes_string(x = 0, y = 0, colour = "target", xend = current.x,
                              yend = current.y),
                   arrow = arrow(length = unit(0.1, "in")), size = 0.1) +
      scale_color_manual(values = c('firebrick2','chartreuse3','dodgerblue4')) 
      # geom_text(data = arrows, aes_string(x = current.x, y = current.y,
      #                                     label = "experiments"), size = 2)

    print(g)
  }
  dev.off()
}
