library(UniProt.ws)

set.name.of.group <- function(protein.groups.f, groups.list, subgroup.list, id.column = 'Protein.IDs') {
  protein.groups.f$group <- "no group"
  protein.groups.f$subgroup <- "no subgroup"
  no.reverse.pattern <- "(?<!(REV__))"
  for (group in names(groups.list)) {
    rows <-
      grep(paste0(no.reverse.pattern,
                  paste(groups.list[[group]], collapse = paste0('|', no.reverse.pattern))),
           unlist(protein.groups.f[, id.column, with = F]), perl = T)
    if (!identical(protein.groups.f[rows, ]$group, character(0))) {
      protein.groups.f[rows, ]$group <- group
    }
  }
  for (subgroup in names(subgroup.list)) {
    rows <-
      grep(paste0(no.reverse.pattern,
                  paste(subgroup.list[[subgroup]], collapse = paste0('|', no.reverse.pattern))),
           unlist(protein.groups.f[, id.column, with = F]), perl = T)
    if (!identical(protein.groups.f[rows, ]$subgroup, character(0))) {
      protein.groups.f[rows, ]$subgroup <- subgroup
    }
  }
  return(protein.groups.f)
}

make.uniprot.gene.names <- function(metadata) {
  if (file.exists(kPathToUniprotGeneNames)) {
    uniprot.gene.names <- 
      fread(kPathToUniprotGeneNames, header = TRUE, sep = "\t", verbose = TRUE)
    row.order <- match(metadata$Protein.IDs, uniprot.gene.names$Protein.IDs)
    metadata$uniprot.genes <- uniprot.gene.names[row.order, uniprot.genes]
  } else {
    up <- UniProt.ws(taxId=9606)
    most.likely.proteins <- unlist(lapply(metadata$Protein.IDs, function(protein.str) {
      strsplit(protein.str, ";", fixed = T)[[1]][1]
    }))
    
    most.likely.proteins <- str_remove(most.likely.proteins, "-[0-9]")
    metadata$uniprot.genes <-
      UniProt.ws::select(up, keys = most.likely.proteins, columns = c("GENES"),
                         keytype = "UNIPROTKB")[[2]]
    metadata$uniprot.genes <-
      unlist(lapply(strsplit(metadata$uniprot.genes, " ", fixed = T), function(str) {
        paste(str, collapse = ";")
      }))
    
    metadata[uniprot.genes == "NA", ]$uniprot.genes <- NA
    metadata[is.na(uniprot.genes), ]$uniprot.genes <- metadata[is.na(uniprot.genes), ]$Gene.names
    
    write.table(metadata[, c("Protein.IDs", "uniprot.genes")], 
                file = kPathToUniprotGeneNames, quote = F, sep = "\t",
                row.names = F)
  }
  
  metadata$uniprot.gene.best <- unlist(lapply(metadata$uniprot.genes, function(genes.str) {
    strsplit(genes.str, ";", fixed = T)[[1]][1]
  }))

  metadata
}

make.metadata.proteins <- function(path.to.protein.groups, groups.list, subgroup.list) {
  prot.groups <- fread(path.to.protein.groups, header = TRUE, sep = "\t", verbose = TRUE)
  
  # de-contaminant GFP
  prot.groups[grep(x = prot.groups$`Protein IDs`, pattern = "Q9U6Y5")]$`Potential contaminant` <- ""
  
  metadata <- prot.groups[, c(1:12), with = F]
  metadata[, potential.contaminant := prot.groups$`Potential contaminant`]
  metadata[, reverse := prot.groups$Reverse]
  metadata[, contaminant.bool := ifelse(potential.contaminant == "+", T, F)]
  metadata[, reverse.bool := ifelse(reverse == "+", T, F)]
  colnames(metadata) <- make.names(colnames(metadata))
  # fill metadata with protein groups
  metadata <-
    set.name.of.group(protein.groups.f = metadata, groups.list = groups.list,
                      subgroup.list = subgroup.list,
                      id.column = 'Protein.IDs')
  metadata <- make.uniprot.gene.names(metadata)
  
  return(metadata)
}


make.metadata.experiments <- function(experiment.names) {
  locations <- vapply(experiment.names, function(current.string) {
    str_remove(current.string, "_.*")
  }, 
  FUN.VALUE = "", USE.NAMES = F
  )
  
  targets <- vapply(experiment.names, function(current.string) {
    str_remove(str_remove(current.string, "[^_]+_"), "_.*")
  }, 
  FUN.VALUE = "", USE.NAMES = F
  )
  
  is.case <- vapply(experiment.names, function(current.string) {
    if (str_remove(str_remove(current.string, "[^_]+_[^_]+_"), "_.*") == "C") {
      return(F)
    }
    
    return(T)
  }, 
  FUN.VALUE = T, USE.NAMES = F
  )
  
  buffers <- vapply(experiment.names, function(current.string) {
    str_remove(str_remove(current.string, ".*_"), "[a-z]")
  }, 
  FUN.VALUE = "", USE.NAMES = F
  )
  
  replicates <- vapply(experiment.names, function(current.string) {
    str_remove(current.string, ".*_[0-9]+")
  }, 
  FUN.VALUE = "", USE.NAMES = F
  )
  
  replicates[grep(pattern = "_GA", x = replicates)] <- 
    str_remove(replicates[grep(pattern = "_GA", x = replicates)], "_GA")
  
  buffers[buffers == "GA"] <- "12_GA"
  
  return(data.table(experiment.name = experiment.names, location = locations, target = targets,
                    is.case = is.case, buffer = buffers, replicate = replicates))
}
