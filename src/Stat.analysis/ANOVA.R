library(data.table)
library(stringr)
library(bit64)
library(ggplot2)
library(ggrepel)
source("./lib/constants.R")

MaxQuant$set("public", "PerformANOVA", function(is.log.norm = T, save.result.table = F, 
                                                do.plot = F, color.list = kColorList, prefix = '',
                                                plot.pvalues = F) {
  private$metadata.proteins$GENES <- private$metadata.proteins$uniprot.genes
  private$metadata.proteins[private$metadata.proteins$uniprot.genes == "" |
                              is.na(private$metadata.proteins$uniprot.genes)]$GENES <-
    private$metadata.proteins[private$metadata.proteins$uniprot.genes == "" |
                                is.na(private$metadata.proteins$uniprot.genes)]$Protein.IDs
  private$metadata.proteins$GENES <- 
    str_split_fixed(private$metadata.proteins$GENES, ';', Inf)[, 1]
  private$metadata.proteins <- 
    private$metadata.proteins[private$metadata.proteins$Razor...unique.peptides > 1, ]
  private$intensity.table <- private$intensity.table[private$intensity.table$`Protein IDs` %in% 
                                                       private$metadata.proteins$Protein.IDs,]
  results <- data.table("Protein.IDs" = private$metadata.proteins$Protein.IDs, 
                        "GENES" = private$metadata.proteins$GENES, 
                        "log2.fold.change" = 0, "pvalue.adj" = - 1)
  n.significant <- list()
  n.known.significant <- list()
  date <- MakeDate()
  if (is.null(private$metadata.experiments$group)){
    private$metadata.experiments$group <- paste(private$metadata.experiments$location, 
                                                private$metadata.experiments$target, 
                                                private$metadata.experiments$buffer, sep = '_')
  } 
  
  raw.pvalues <- list()
  
  if (do.plot)
    pdf(paste('out/', date, "volcano_plots.pdf", sep = ""))
  for(exp in unique(private$metadata.experiments$group)){
    print(exp)
    a <- private$metadata.experiments[private$metadata.experiments$group == exp &
                                      private$metadata.experiments$is.case == T, ]$experiment.name
    b <- private$metadata.experiments[private$metadata.experiments$group == exp &
                                      private$metadata.experiments$is.case == F, ]$experiment.name
    pr.f <- private$intensity.table[, c("Protein IDs", a, b), with = F]
    pr.f <- merge(private$metadata.proteins[, c('Protein.IDs', 'GENES', 'group')], pr.f, 
                  by.x = 'Protein.IDs', by.y = 'Protein IDs')
    pr.f <- pr.f[rowSums(!is.na(pr.f)) > 3] # remove not identified proteins
    pr.f[["pvalue"]] <- NA
    for (i in 1:nrow(pr.f)) {
      if (length(unique(unlist(pr.f[i, a, with = F]))) == 1 | 
          length(unique(unlist(pr.f[i, b, with = F]))) == 1) {
        pr.f[["pvalue"]][i] <- 0  
      } else {
        pr.f[["pvalue"]][i] <- 
          t.test(unlist(pr.f[i, a, with = F]), 
                 unlist(pr.f[i, b, with = F]), var.equal = T)$p.value  
      }
      
      mean.a <- mean(unlist(pr.f[i, a, with = F]))
      mean.b <- mean(unlist(pr.f[i, b, with = F]))
      if (is.log.norm) {
        pr.f[["log2.fold.change"]][i] <- mean.a - mean.b  
      } else {
        pr.f[["log2.fold.change"]][i] <- mean.a/mean.b  
      }
      
    }
    
    pr.f[["pvalue.adj"]] <- p.adjust(pr.f[["pvalue"]], method = "BH")
    drops <- c(a, b)
    pr.f <- pr.f[, -drops, with = F]
    n.significant[[exp]] <- nrow(pr.f[pr.f$pvalue.adj < 0.05 & pr.f$log2.fold.change > 1, ])
    n.known.significant[[exp]] <- nrow(pr.f[pr.f$pvalue.adj < 0.05 & pr.f$log2.fold.change > 1 & 
                                              pr.f$group != 'no group',])
    # volcano plot
    if (do.plot)

      print(ggplot(as.data.frame(pr.f), aes(x = log2.fold.change, y = -log10(pvalue.adj),
                                            color = group)) +
              geom_rect(aes(xmin = 1, xmax = Inf, ymin = -log10(.05), ymax = Inf),
                        fill = "gray85", alpha = 0.03, linetype = "blank") +
              geom_point(aes(color = group), size = 3) +
              scale_color_manual(values = color.list) +
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
              ggtitle(exp) +
              theme(legend.direction='vertical')) 
    
    raw.pvalues[[exp]] <- pr.f$pvalue
    # save
    drops <- c("group", 'pvalue')
    pr.f <- pr.f[, -drops, with = F]
    results <- merge(results, pr.f, all = T, by = c("Protein.IDs", "GENES"), 
                     suffixes = c("",paste(".",exp,sep = '')))
  }
  
  if (do.plot)
    dev.off()
  drops <- c("log2.fold.change", "pvalue.adj")
  results <- results[, !(names(results) %in% drops), with = F]
  if (save.result.table) {
    saveRDS(results, paste0('out/', prefix, MakeDate(), '_pvals.logfc.rds'))
  }
  
  if (plot.pvalues) {
    pdf(paste0('out/', MakeDate(), "_pvalues_distribution.pdf"), width = 16, height = 16)
    
    layout(matrix(c(1:length(raw.pvalues)), ncol = 4, byrow = F))
    
    for (exp in names(raw.pvalues)) {
      hist(raw.pvalues[[exp]], main = exp, col = "black", xlab = "raw p-values",
           breaks = 100)
    }
    
    dev.off()
  }
  
  private$pvals.logfc <- results
  private$AddChange("performed ANOVA and calculated logFC")
  return(list(n.significant, n.known.significant))  
  invisible(self)
})

