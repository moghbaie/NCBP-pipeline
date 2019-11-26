require(data.table)
require(stringr)
require(R6)

source("lib/constants.R")
source("src/Preprocessing/rename.files.R")
source("src/Preprocessing/make.metadata.R")

MaxQuant <- R6Class("MaxQuant",
                    private = list(
                      intensity.table = NULL,
                      metadata.proteins = NULL,
                      metadata.experiments = NULL,
                      pvals.logfc = NULL,
                      distance.matrix = NULL,
                      clusters = NULL,
                      changes = c(), 
                      changes.count = 0,
                      is.zero.intensity = NULL,
                      hypotheses.status = NULL,
                      
                      AddChange = function(change.string) {
                        if (is.null(change.string)) {
                          stop("provide a change string")
                        }
                        private$changes.count = private$changes.count + 1
                        private$changes[[private$changes.count]] = change.string
                      }
                    ),
                    public = list(
                      initialize = function(intensity.table, metadata.proteins, 
                                            metadata.experiments, pvals.logfc = NULL,
                                            changes = NULL) {
                        private$intensity.table <- intensity.table
                        private$metadata.proteins <- metadata.proteins
                        private$metadata.experiments <- metadata.experiments
                        private$pvals.logfc <- pvals.logfc
                        if (length(changes) > 0) {
                          private$changes <- changes
                          private$changes.count <- length(changes)
                        }
                      },
                      print = function() {
                        cat("List of changes:\n")
                        i = 1
                        for (change.string in private$changes) {
                          cat(paste0(i, "."), change.string, sep = " ")
                          i <- i + 1
                          cat("\n")
                        }
                      },
                      SetIntensityTable = function(new.intensity.table,
                                                     change.string = NULL) {
                        private$AddChange(change.string)
                        private$intensity.table <- new.intensity.table
                      },
                      GetIntensityTable = function() {
                        return(private$intensity.table)
                      }, 
                      SetMetadataProteins = function(new.metadata.proteins,
                                                       change.string = NULL) {
                        private$AddChange(change.string)
                        private$metadata.proteins <- new.metadata.proteins
                      },
                      GetMetadataProteins = function() {
                        return(private$metadata.proteins)
                      }, 
                      SetMetadataExperiments = function(new.metadata.experiments,
                                                          change.string = NULL) {
                        private$AddChange(change.string)
                        private$metadata.experiments <- new.metadata.experiments
                      },
                      GetMetadataExperiments = function() {
                        return(private$metadata.experiments)
                      },
                      GetPvaluesLogFC = function() {
                        return(private$pvals.logfc)
                      },
                      GetDistanceMatrix = function() {
                        return(private$distance.matrix)
                      },
                      GetClusters = function() {
                        return(private$clusters)
                      },
                      GetIsZeroIntensity = function() {
                        return(private$is.zero.intensity)
                      },
                      GetHypothesesStatus = function() {
                        return(private$hypotheses.status)
                      }
                     ))

BuildMaxQuant <- function(path.to.protein.groups, pattern = "iBAQ",
                          groups.list, subgroup.list, new.data = T) {
  prot.groups <- fread(path.to.protein.groups, header = TRUE, sep = "\t", verbose = TRUE)
  
  prot.groups <- prot.groups[, grep(x = colnames(prot.groups), pattern = paste0(pattern, ".+")), 
                             with = F]
  prot.groups <- rename.files.to.experiments(protein.groups.f = prot.groups, 
                                             pattern. = pattern, new.data = new.data)
  metadata.proteins <- make.metadata.proteins(path.to.protein.groups = path.to.protein.groups, 
                                              groups.list = groups.list,
                                              subgroup.list = subgroup.list)
  metadata.experiments <- make.metadata.experiments(colnames(prot.groups))
  prot.groups[[kIdColumnName]] <- metadata.proteins$Protein.IDs
  return(MaxQuant$new(prot.groups, metadata.proteins, metadata.experiments))
}

