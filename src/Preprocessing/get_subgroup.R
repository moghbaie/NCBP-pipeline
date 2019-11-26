require(data.table)
source("lib/helper_functions.R")

MaxQuant$set("public", "MakeObjectByGroupName", function(group_name) {
  if (is.null(private$metadata.proteins)) {
    stop("You don't have metadata file")
  }
  
  groups <- unique(private$metadata.proteins$group)
  if (group_name %nin% groups) {
    stop("Selected group name out of range")
  }
  
  group.protein.ids <- private$metadata.proteins[group == group_name]$Protein.IDs
  instensities.subgroup <- private$intensity.table[`Protein IDs` %in% group.protein.ids]
  if (!is.null(private$pvals.logfc)) {
    pvals.logfc.subgroup <- private$pvals.logfc[`Protein.IDs` %in% group.protein.ids]
  } else {
    pvals.logfc.subgroup <- NULL
  }
  metadata.proteins.subgroup <- private$metadata.proteins[group == group_name]
  
  return(MaxQuant$new(instensities.subgroup, metadata.proteins.subgroup, 
                      private$metadata.experiments, pvals.logfc.subgroup, 
                      c(private$changes, paste0("Subgroup selected: ", group_name))))
})
