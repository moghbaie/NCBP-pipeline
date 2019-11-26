library(grid)
library(ggpubr)


PlotVolcanoPlots <- function(max.quant.obj, mode = 'all'){
  dt <- max.quant.obj$GetPvaluesLogFC()
  metadata.experiment <- max.quant.obj$GetMetadataExperiments()
  metadata.experiment$group <- paste(metadata.experiment$location, metadata.experiment$target, metadata.experiment$buffer, sep = '_')
  metadata.protein <- max.quant.obj$GetMetadataProteins()
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Gene.names
  metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$uniprot.gene.best <- metadata.protein[is.na(metadata.protein$uniprot.gene.best)]$Protein.IDs
  metadata.protein$uniprot.gene.best <- str_split_fixed(metadata.protein$uniprot.gene.best, ';', 2)[, 1]
  dt <- merge(dt, metadata.protein[, c('Protein.IDs', 'group', 'uniprot.gene.best'), with = F], by = 'Protein.IDs', all.x = T)
  dt$GENES <- dt$uniprot.gene.best
  pdf(paste('out/', MakeDate(), "volcano_plots.pdf", sep = ""), width = 10, height = 10)
  
  for(exp in unique(metadata.experiment$group)){
    print(exp)
    pretty_name <- paste(unique(metadata.experiment[metadata.experiment$group == exp, ]$pretty_name), 
                         unique(metadata.experiment[metadata.experiment$group == exp, ]$buffer), sep = '_')
    
    pr.f <- dt[, c('Protein.IDs', 'GENES', 'group', grep(exp, colnames(dt), value = T)), with = F]
    colnames(pr.f) <- c('Protein.IDs', 'GENES', 'group', 'log2.fold.change', 'pvalue.adj')
    if(mode == 'all'){
      g <- ggplot(pr.f, aes(x = log2.fold.change, y = -log10(pvalue.adj),
                            color = group)) +
        geom_rect(aes(xmin = 1, xmax = Inf, ymin = -log10(.05), ymax = Inf),
                  fill = "gray85", alpha = 0.03, linetype = "blank") +
        geom_point(aes(color = group), size = 3) +
        scale_color_manual(values = kColorList) +
        theme_bw() +
        geom_text_repel(
          data = pr.f[pr.f$group != 'no group', ],
          aes(label = GENES),
          size = 3,
          box.padding = unit(0.4, "lines"),
          point.padding = unit(0.4, "lines"))   +
        geom_hline(yintercept = -log10(.05), color = "grey", linetype="dashed") +
        geom_vline(xintercept = -1, color = "grey", linetype="dashed") +
        geom_vline(xintercept = 1, color = "grey", linetype="dashed") +
        ggtitle(pretty_name) +
        theme(legend.direction='vertical')   
    } else if(mode == 'no_spliceosome'){
      pr.f[(pr.f$log2.fold.change <= 1 | pr.f$pvalue.adj >= 0.05) & !grepl(paste(kTargetProteins, collapse = '|'), pr.f$Protein.IDs), ]$group <- 'no group'
      g <- ggplot(pr.f, aes(x = log2.fold.change, y = -log10(pvalue.adj),
                            color = group)) +
        geom_rect(aes(xmin = 1, xmax = Inf, ymin = -log10(.05), ymax = Inf),
                  fill = "gray85", alpha = 0.03, linetype = "blank") +
        geom_point(aes(color = group), size = 3) +
        scale_color_manual(values = kColorList) +
        theme_bw() +
        geom_text_repel(
          data = pr.f[(pr.f$group != 'no group' & 
                       pr.f$log2.fold.change > 1 & 
                       pr.f$pvalue.adj < 0.05 & 
                       pr.f$group != "Spliceosome") | 
                      ((pr.f$group == 'NCBP1/2-related' | 
                        pr.f$group == 'NCBP3')), ],
          aes(label = GENES),
          size = 3,
          box.padding = unit(0.4, "lines"),
          point.padding = unit(1, "lines"))   +
        geom_hline(yintercept = -log10(.05), color = "grey", linetype="dashed") +
        geom_vline(xintercept = -1, color = "grey", linetype="dashed") +
        geom_vline(xintercept = 1, color = "grey", linetype="dashed") +
        ggtitle(pretty_name) +
        theme(legend.direction='vertical')   
    } else if(mode == 'outside'){
      texts.list <- list()
      texts.list[1] <- list(NULL) 
      counts.list <- list()
      for(group. in unique(pr.f[pr.f$group != 'no group', ]$group)){
        text <- paste(pr.f[pr.f$log2.fold.change > 1 & pr.f$pvalue.adj < 0.05 & pr.f$group == group., ]$GENES, collapse = ', ')
        if(text != '' & !grepl(',', text)){
          text <- paste0(' ', text, ' ')
        }
        if(text != ''){
          text.p <- ggparagraph(text = text, size = 7, color = kColorList[group.])
          texts.list[[group.]] <- text.p  
          counts.list[[group.]] <- length(pr.f[pr.f$log2.fold.change > 1 & pr.f$pvalue.adj < 0.05 & pr.f$group == group., ]$GENES)
        }
      }
      l <- length(texts.list)
      texts.list[l + 1] <- list(NULL) 
      h <- c(1, 0.1 * (unname(unlist(counts.list)) %/% 4 + 1), 1)
      texts <- ggarrange(plotlist = texts.list, 
                         ncol = 1, nrow = length(texts.list),
                         heights = h)
      pr.f[(pr.f$log2.fold.change <= 1 | pr.f$pvalue.adj >= 0.05) & !grepl(paste(kTargetProteins, collapse = '|'), pr.f$Protein.IDs), ]$group <- 'no group'
      g <- ggplot(pr.f, aes(x = log2.fold.change, y = -log10(pvalue.adj),
                            color = group)) +
        geom_rect(aes(xmin = 1, xmax = Inf, ymin = -log10(.05), ymax = Inf),
                  fill = "gray85", alpha = 0.03, linetype = "blank") +
        geom_point(aes(color = group), size = 3) +
        scale_color_manual(values = kColorList) +
        theme_bw() +
        geom_hline(yintercept = -log10(.05), color = "grey", linetype="dashed") +
        geom_vline(xintercept = -1, color = "grey", linetype="dashed") +
        geom_vline(xintercept = 1, color = "grey", linetype="dashed") +
        ggtitle(pretty_name) +
        theme(legend.direction='horizontal',
              legend.position = 'top') 
      
      g <- ggarrange(g, texts, 
                ncol = 2, nrow = 1, widths = c(1, 0.3))
      
    }
    print(g)
  }
  dev.off()
}




