library(data.table)
library(stringr)
library(bit64)
library(dplyr)


MaxQuant$set("public", "ControlFPFN", function(passed.anova = F, save.table = T,
                                               save.fp.0.1 = F, mode = "single") {
  
  metadata.experiments <- private$metadata.experiments
  row.order <- match(private$pvals.logfc$Protein.IDs, private$is.zero.intensity[[kIdColumnName]])
  # row.order <- row.order[!is.na(row.order)]
  private$is.zero.intensity <- private$is.zero.intensity[row.order, ]
  
  experimental.groups <- unique(private$metadata.experiments$group)
  
  logfc.cols <- colnames(private$pvals.logfc)[grep("log", colnames(private$pvals.logfc))]
  pval.cols <- colnames(private$pvals.logfc)[grep("pval", colnames(private$pvals.logfc))]
    
  hypotheses.status <- data.table(matrix("N", ncol = length(experimental.groups), 
                                    nrow = dim(private$is.zero.intensity)[[1]]))
  colnames(hypotheses.status) <- experimental.groups
  hypotheses.status$Protein.ID <- private$is.zero.intensity[[kIdColumnName]]
  
  fp_0_1 = list()
    
  for (group.name in experimental.groups) {
    case.cols <- private$metadata.experiments[(group == group.name) & is.case, experiment.name]
    control.cols <-
      private$metadata.experiments[(group == group.name) & (!is.case), experiment.name]
      
    # check if protein is significant
    logfc.col <- logfc.cols[grep(paste0(group.name, "$"), logfc.cols)]
    pval.col <- pval.cols[grep(paste0(group.name, "$"), pval.cols)]
      
    if (mode == "single") {
      is.signif.experiment <-
        ((private$pvals.logfc[, ..pval.col] < 0.05) & (private$pvals.logfc[, ..logfc.col] > 1))
      # if yes --- FP-1, FP-0 or TP?
      is.signif.experiment[is.na(is.signif.experiment)] <- F
    } else if (mode == "vote") {
      significant.dt <- fread('out/significant_proteins.tsv')
      proteins.sign <- significant.dt[experiment == group.name, Protein.IDs]
      rm(significant.dt)
      
      is.signif.experiment <- unlist(lapply(private$pvals.logfc$Protein.IDs, function(protein) {
        protein %in% proteins.sign
      }))
    } else {
      stop(paste("Unknown value for mode:", mode))
    }
      
    candidates.fp.1 <- rowSums(private$is.zero.intensity[, case.cols, with = F]) == 1
    candidates.fp.0 <- (rowSums(private$is.zero.intensity[, case.cols, with = F]) == 0) &
      (rowSums(private$is.zero.intensity[, control.cols, with = F]) > 0)
    candidates.fn <- (rowSums(private$is.zero.intensity[, case.cols, with = F]) > 1) &
      (rowSums(private$is.zero.intensity[, control.cols, with = F]) <= 1)
    
    fp1 <- is.signif.experiment & candidates.fp.1 
    fp0 <- is.signif.experiment & candidates.fp.0
    
    fp <- fp1 | fp0
    data.table::set(x = hypotheses.status, i = which(fp), j = group.name, value = "FP")
    
    fn <- !is.signif.experiment & candidates.fn
    data.table::set(x = hypotheses.status, i = which(fn), j = group.name, value = "FN")
    
    tp <- is.signif.experiment & !(candidates.fp.1 | candidates.fp.0)
    data.table::set(x = hypotheses.status, i = which(tp), j = group.name, value = "TP")
      
    tn <- !is.signif.experiment & !candidates.fn
    data.table::set(x = hypotheses.status, i = which(tn), j = group.name, value = "TN")
    
    fp_0_1[[group.name]] = data.table(FP1 = sum(fp1), FP0 = sum(fp0))
  }
  
  fp_0_1  = rbindlist(fp_0_1)
  fp_0_1$experiment = experimental.groups
  
  if (save.fp.0.1) {
    write.table(fp_0_1, file = paste0("out/", MakeDate(), "_fp_0_1.csv"), row.names = F,
                sep = "\t", quote = F)
  }
  
  private$hypotheses.status <- hypotheses.status
  private$AddChange("added hypotheses status")
  
  if (passed.anova) {
    fp.proteins <- unlist(hypotheses.status[(rowSums(hypotheses.status == "FP") != 0) & 
                                              (rowSums(hypotheses.status == "TP") > 0),
                                            "Protein.ID"]) 
  } else {
    fp.proteins <- unlist(hypotheses.status[(rowSums(hypotheses.status == "FP") != 0) & 
                                              (rowSums(hypotheses.status == "TP") == 0),
                                            "Protein.ID"]) 
  }
  
  fp.rows <- match(fp.proteins, private$metadata.proteins$`Protein.IDs`)
  
  dt.fp <- data.table(Protein.ID = fp.proteins, 
                      Gene.name = 
                        unlist(private$metadata.proteins[fp.rows, "Gene.names", with = F]),
                      Group = 
                        unlist(private$metadata.proteins[fp.rows, "group", with = F]))
  
  experiments.fp <- apply(hypotheses.status[Protein.ID %in% fp.proteins, ], 1, function(row) {
    paste(colnames(hypotheses.status)[row == "FP"], collapse = "; ")
  })
  
  dt.fp$`Experiments where FP` = 
    experiments.fp[match(dt.fp$Protein.ID, 
                         hypotheses.status[Protein.ID %in% fp.proteins, Protein.ID])]
  
  experiments.tp <- apply(hypotheses.status[Protein.ID %in% fp.proteins, ], 1, function(row) {
    paste(colnames(hypotheses.status)[row == "TP"], collapse = "; ")
  })
  
  dt.fp$`Experiments where TP` = 
    experiments.tp[match(dt.fp$Protein.ID, 
                         hypotheses.status[Protein.ID %in% fp.proteins, Protein.ID])]
  
  group.order <- c(unique(unlist(kGroupOrder[, "Function"])), kNoGroupName)
  
  dt.fp = dt.fp %>%
    arrange(match(Group, group.order))
  
  if (passed.anova) {
    passed.anova.str = "_passed_anova_once"
  } else {
    passed.anova.str = "_didnt_pass_anova"
  }
  
  dt.fp$`Experiments where FP` <- unlist(lapply(dt.fp$`Experiments where FP`, function(x) {
    exps <- str_split(x, '; ')
    paste(paste(metadata.experiments[match(exps[[1]], metadata.experiments$group), ]$pretty_name, str_split_fixed(exps[[1]], '_', 2)[, 2], sep = '_'), collapse = '; ')
  }))
  
  dt.fp$`Experiments where TP` <- unlist(lapply(dt.fp$`Experiments where TP`, function(x) {
    exps <- str_split(x, '; ')
    paste(paste(metadata.experiments[match(exps[[1]], metadata.experiments$group), ]$pretty_name, str_split_fixed(exps[[1]], '_', 2)[, 2], sep = '_'), collapse = '; ')
  }))
  
  dt.fp$`Experiments where TP` <- gsub('NA_', " ", dt.fp$`Experiments where TP`)
  dt.fp$`Experiments where FP` <- gsub('NA_', " ", dt.fp$`Experiments where FP`)
  
  if (save.table) {
    write.table(dt.fp, file = paste0("out/", MakeDate(), passed.anova.str, "_FP_proteins.csv"), 
                quote = F, sep = "\t", row.names = F)    
  }

  
  if (passed.anova) {
    fn.proteins <- unlist(hypotheses.status[(rowSums(hypotheses.status == "FN") != 0) & 
                                              (rowSums(hypotheses.status == "TP") > 0),
                                            "Protein.ID"]) 
  } else {
    fn.proteins <- unlist(hypotheses.status[(rowSums(hypotheses.status == "FN") != 0) & 
                                              (rowSums(hypotheses.status == "TP") == 0),
                                            "Protein.ID"]) 
  }
  
  fn.rows <- match(fn.proteins, private$metadata.proteins$`Protein.IDs`)
  
  dt.fn <- data.table(Protein.ID = fn.proteins, 
                      Gene.name = 
                        unlist(private$metadata.proteins[fn.rows, "Gene.names", with = F]),
                      Group = 
                        unlist(private$metadata.proteins[fn.rows, "group", with = F]))
  
  experiments.fn <- apply(hypotheses.status[Protein.ID %in% fn.proteins, ], 1, function(row) {
    paste(colnames(hypotheses.status)[row == "FN"], collapse = "; ")
  })
  
  dt.fn$`Experiments where FN` = 
    experiments.fn[match(dt.fn$Protein.ID, 
                         hypotheses.status[Protein.ID %in% fn.proteins, Protein.ID])]
  
  experiments.tp <- apply(hypotheses.status[Protein.ID %in% fn.proteins, ], 1, function(row) {
    paste(colnames(hypotheses.status)[row == "TP"], collapse = "; ")
  })
  
  dt.fn$`Experiments where TP` = 
    experiments.tp[match(dt.fn$Protein.ID, 
                         hypotheses.status[Protein.ID %in% fn.proteins, Protein.ID])]
  
  dt.fn = dt.fn %>%
    arrange(match(Group, group.order))

  
  dt.fn$`Experiments where FN` <- unlist(lapply(dt.fn$`Experiments where FN`, function(x) {
    exps <- str_split(x, '; ')
    paste(paste(metadata.experiments[match(exps[[1]], metadata.experiments$group), ]$pretty_name, str_split_fixed(exps[[1]], '_', 2)[, 2], sep = '_'), collapse = '; ')
  }))
  
  dt.fn$`Experiments where TP` <- unlist(lapply(dt.fn$`Experiments where TP`, function(x) {
    exps <- str_split(x, '; ')
    paste(paste(metadata.experiments[match(exps[[1]], metadata.experiments$group), ]$pretty_name, str_split_fixed(exps[[1]], '_', 2)[, 2], sep = '_'), collapse = '; ')
  }))
  
  dt.fn$`Experiments where TP` <- gsub('NA_', " ", dt.fn$`Experiments where TP`)
  dt.fn$`Experiments where FN` <- gsub('NA_', " ", dt.fn$`Experiments where FN`)
  
  if (save.table) {
    write.table(dt.fn, file = paste0("out/", MakeDate(), passed.anova.str, "_FN_proteins.csv"), 
                quote = F, sep = "\t", row.names = F)    
  }
  
  return(list("#FP" = length(fp.proteins), "#FN" = length(fn.proteins)))
})



