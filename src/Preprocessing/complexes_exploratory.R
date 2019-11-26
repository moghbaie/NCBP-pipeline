source("lib/helper_functions.R")

require(data.table)
require(stringr)

path.to.complexes <- "data/complexes.csv"
complexes <- fread(path.to.complexes, header = TRUE, sep = "\t", verbose = TRUE)

# merge spliceosome

spliceosome.name <- "Spliceosome"

spliceosome.proteins <-
  unique(strsplit(paste(complexes[Function == spliceosome.name, ]$uniprotID, collapse = "|"), 
                "|", fixed = T)[[1]])

complexes <- complexes[Function != spliceosome.name, ]
complexes <-
  rbindlist(list(complexes, data.table(Function = spliceosome.name, ComplexName = spliceosome.name,
                                     uniprotID = paste(spliceosome.proteins, collapse = "|"))))

# check for intersection

complexes.long <- 
  complexes[, .(uniprotID = strsplit(uniprotID, "|", fixed = T)[[1]]), by = c("Function", 
                                                                            "ComplexName")]
complexes.long <- unique(complexes.long, by = c("uniprotID"))

write.table(complexes.long, file = paste0("data/corum_complexes_", MakeDate(), ".csv"), 
            quote = F, sep = "\t", row.names = F)
