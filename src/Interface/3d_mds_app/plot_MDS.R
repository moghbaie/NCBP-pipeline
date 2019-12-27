require(plotly)
require(ape)
require(vegan)
require(ggplot2)

#Sys.setenv("plotly_username"="KalmSveta")
#Sys.setenv("plotly_api_key"="L4QGJcqZsrNTVySXjs3W")

plotly.color.list <- c("grey", "red", "dark green", "orange", "turquoise4", "green", "blue", 
                       "black", "brown", "violet", "cyan")
names(plotly.color.list) <- c("no group", "NCBP1/2-related", "EJC-related", "THO/TREX-related", 
                              "Exosome-related", "NPC-related", "ERH", "NELF",
                              "Spliceosome", "Importin related", "NCBP3")
kNoGroupName <- "no group"
'%nin%' <- Negate('%in%')


CountPcoaDf <- function(dist.mtx) {
  pcoa.obj <- pcoa(dist.mtx)
  pcoa.df <- as.data.frame(pcoa.obj$vectors)
  pcoa.df <- pcoa.df[,c(1,2,3)]
  pcoa.df$Protein.IDs <- rownames(pcoa.df)
  pcoa.df
}

Plot3DMDS <- function(max.quant.obj, pcoa.df = NULL, title = "MDS",
                      color.list = plotly.color.list,
                      anovas.passed = 1, groups.to.show = c(), recalculate.dist = F) {
  
  
  metadata.table <- 
    max.quant.obj$GetMetadataProteins()[,c('Protein.IDs','Gene.names','group')]
  metadata.table[metadata.table$group=="NCBP3","Gene.names"] <- "NCBP3"
  metadata.table[metadata.table$group=="Exosome-connected","group"] <- "Exosome-related"
  metadata.table[metadata.table$group=="ERH-related","group"] <- "ERH"
  
  
  df_sign <- max.quant.obj$GetPvaluesLogFC()
  proteins.to.plot <- intersect(
    unlist(df_sign[df_sign$n.of.anovas.passed >= anovas.passed, "Protein.IDs"]),
    unlist(metadata.table[group %in% groups.to.show, 'Protein.IDs']))
  
  validate(need(length(proteins.to.plot) > 3, label = 'More than 3 peptides'))
  
  if (recalculate.dist) {
    # pcoa.df <- CountPcoaDf(subset(max.quant.obj$GetDistanceMatrix(), proteins.to.plot))
    pcoa.df <- 
      CountPcoaDf(as.dist(as.matrix(max.quant.obj$GetDistanceMatrix())[
        proteins.to.plot, proteins.to.plot]))
  } else {
    if (is.null(pcoa.df)) {
      pcoa.df <- CountPcoaDf(max.quant.obj$GetDistanceMatrix())
    }
    
    pcoa.df <- pcoa.df[pcoa.df$Protein.IDs %in% proteins.to.plot, ]
  }
  
  pcoa.df <- merge(pcoa.df, metadata.table, by = 'Protein.IDs', all.x = T, all.y = F)
  pcoa.df$always_names <- ''
  pcoa.df[pcoa.df$group != kNoGroupName, ]$always_names <- 
    pcoa.df[pcoa.df$group != kNoGroupName,]$Gene.names
  
  p <- 
    plot_ly(pcoa.df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, 
            color = ~group, colors = color.list,
            marker = list(size = 2), text = ~always_names) %>%
    add_markers(showlegend = F) %>% 
    add_text(textposition = "top right", showlegend = F) %>%
    layout(scene = list(xaxis = list(title = 'Axis.1'),
                        yaxis = list(title = 'Axis.2'),
                        zaxis = list(title = 'Axis.3'))) %>%
    add_trace(x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
              text = ~Gene.names,
              hoverinfo = 'text'
    ) %>%
    layout(title = title)
}
