require(data.table)

MaxQuant$set("public", "MakeObjectOnlySignificant", function(input.pvals.logfc = NULL, mode = 'vote') {
  if(mode == 'single'){
    if (is.null(private$pvals.logfc) & is.null(input.pvals.logfc)) {
      stop("You need to perform ANOVA first or provide a pvals.logfc table")
    }
    groups <- unique(private$metadata.experiments$group)
    
    if (is.null(input.pvals.logfc)) {
      pvals.logfc <- private$pvals.logfc
    } else {
      pvals.logfc <- input.pvals.logfc
    }
    
    logfc.cols <- 
      colnames(private$pvals.logfc)[grep("log2", colnames(pvals.logfc))]
    pval.cols <- 
      colnames(private$pvals.logfc)[grep("pval", colnames(pvals.logfc))]
    
    list.of.significant.per.column <- lapply(groups, function(group) {
      (pvals.logfc[, paste0("log2.fold.change.", group), with = F] > 1) &
        (pvals.logfc[, paste0("pvalue.adj.", group), with = F] < 0.05)
    })
    
    sign.m <- do.call("cbind", list.of.significant.per.column)
    rownames(sign.m) <- pvals.logfc$Protein.IDs
    colnames(sign.m) <- groups
    sign.m <- replace(sign.m, is.na(sign.m), F)
    
    n.of.anovas.passed <- rowSums(sign.m == T)
    n.of.anovas.passed <- n.of.anovas.passed[n.of.anovas.passed > 0]
    
    row.sign <- match(names(n.of.anovas.passed), private$intensity.table[[kIdColumnName]])
    intensities.sign <- private$intensity.table[row.sign, ]
    intensities.sign[is.na(intensities.sign)] <- 0  
    
    rows.sign.metadata.proteins <- 
      match(names(n.of.anovas.passed), private$metadata.proteins$Protein.IDs)
    metadata.proteins.sign <- private$metadata.proteins[rows.sign.metadata.proteins, ]
    
    rows.sign.pvals.logfc <- match(names(n.of.anovas.passed), pvals.logfc$Protein.IDs)
    pvals.logfc.sign <- pvals.logfc[rows.sign.pvals.logfc, ]
    pvals.logfc.sign$n.of.anovas.passed <- 
      n.of.anovas.passed[match(names(n.of.anovas.passed), pvals.logfc.sign$Protein.IDs)]
  } else if (mode == 'vote') {
    significant.dt <- fread('out/significant_proteins.tsv')
    proteins.sign <- significant.dt$Protein.IDs
    intensities.sign <- private$intensity.table[private$intensity.table$`Protein IDs` %in% proteins.sign, ]
    intensities.sign[is.na(intensities.sign)] <- 0  
    metadata.proteins.sign <- private$metadata.proteins[private$metadata.proteins$Protein.IDs %in% proteins.sign, ]
    n.of.anovas.passed <- as.data.frame(table(significant.dt$Protein.IDs))
    pvals.logfc <- fread('out/pvalue_logfc.tsv')
    pvals.logfc <- pvals.logfc[, names(pvals.logfc) %nin% c('sd'), with = F]
    pvals.logfc <- dcast(pvals.logfc, Protein.IDs + GENES ~ stat + experiment, value.var = 'mean')
    pvals.logfc.sign <- pvals.logfc[pvals.logfc$Protein.IDs %in% proteins.sign, ]
    pvals.logfc.sign <- merge(pvals.logfc.sign, n.of.anovas.passed, by.x = 'Protein.IDs', by.y = 'Var1', all.x = T)
    colnames(pvals.logfc.sign)[colnames(pvals.logfc.sign) == 'Freq'] <- 'n.of.anovas.passed'
  }
  return(MaxQuant$new(intensities.sign, metadata.proteins.sign, 
                      private$metadata.experiments, pvals.logfc.sign, 
                      c(private$changes, 
                        "dropped all proteins, that didn't pass ANOVA at least once\n
                        added column n.of.anovas.passed to the intensity table")))
})
