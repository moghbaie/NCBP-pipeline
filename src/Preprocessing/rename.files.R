rename.files.to.experiments = function(protein.groups.f, pattern. = "LFQ intensity",
                                       new.data = T) {
  
  if (!new.data) {
    lfq.columns <- grep(x = colnames(protein.groups.f), pattern = pattern., value = T)
    setnames(protein.groups.f,
             lfq.columns,
             gsub(".* (.*)", "\\1",
                  lfq.columns))
    colnames(protein.groups.f) <-
      gsub("NCBP2_3", "NCBP23", colnames(protein.groups.f))
    colnames(protein.groups.f) <-
      gsub("tk_", "Columbia_", colnames(protein.groups.f))
    colnames(protein.groups.f) <-
      gsub("Cl_", "C_", colnames(protein.groups.f))
    colnames(protein.groups.f) <-
      gsub("Control_", "NCBP3_C_", colnames(protein.groups.f))
    x = colnames(protein.groups.f[, grep("Columbia_NCBP3_[0-9]", 
                                         colnames(protein.groups.f)), with = F])
    colnames(protein.groups.f)[grep("Columbia_NCBP3_[0-9]", colnames(protein.groups.f))] <-
      unlist(lapply(x, function(x) {
        paste(
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[1],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[2],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[4],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[3],
          sep = "_"
        )
      }))
    x = colnames(protein.groups.f[, grep("Columbia_NCBP3_C", colnames(protein.groups.f)), with = F])
    colnames(protein.groups.f)[grep("Columbia_NCBP3_C", colnames(protein.groups.f))] <-
      unlist(lapply(x, function(x) {
        paste(
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[1],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[2],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[3],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[5],
          unlist(strsplit(
            x, split = "_", fixed = T
          ))[4],
          sep = "_"
        )
      }))
    exps =  colnames(protein.groups.f[, grep("Columbia_NCBP3_[0-9]", 
                                             colnames(protein.groups.f)), with = F])
    colnames(protein.groups.f)[grep("Columbia_NCBP3_[0-9]", colnames(protein.groups.f))] <-
      unlist(lapply(exps, function(exp)
        gsub(
          "(Columbia_NCBP3_.*)(_[0-9])",
          paste("\\1", letters[as.numeric(unlist(strsplit(exp, "_", fixed =
                                                            T))[4])], sep = ""),
          exp
        )))
    exps =  colnames(protein.groups.f[, grep("Columbia_NCBP3_C", colnames(protein.groups.f)), 
                                      with = F])
    colnames(protein.groups.f)[grep("Columbia_NCBP3_C", colnames(protein.groups.f))] <-
      unlist(lapply(exps, function(exp)
        gsub(
          "(Columbia_NCBP3_C_.*)(_[0-9])",
          paste("\\1", letters[as.numeric(unlist(strsplit(exp, "_", fixed = T))[5])], sep = ""),
          exp
        )))
    colnames(protein.groups.f) <-
      gsub("^NCBP", "Rockefeller_NCBP", colnames(protein.groups.f))
    colnames(protein.groups.f)[grep("NCBP23", colnames(protein.groups.f))] <-
      gsub("NCBP23", "NCBP2", colnames(protein.groups.f)[grep("NCBP23", 
                                                             colnames(protein.groups.f))])
    
    df <-
      protein.groups.f[, grep("NCBP2_C", colnames(protein.groups.f)), with = F]
    colnames(df) <- gsub("NCBP2", "NCBP3", colnames(df))
    if (length(df) != 0) protein.groups.f <- cbind(protein.groups.f, df)
    if (any(colnames(protein.groups.f) == pattern.)) protein.groups.f[, (pattern.) := NULL]
  } else {
    colnames(protein.groups.f) <- 
      str_remove(colnames(protein.groups.f), pattern = paste0(pattern., " "))
    
    NCBP12.cols <- grep(pattern = "NCBP[1|2]", x = colnames(protein.groups.f))
    colnames(protein.groups.f)[NCBP12.cols] <- 
      paste0("Rockefeller_", colnames(protein.groups.f)[NCBP12.cols])
    
    NCBP3.cols <- grep(pattern = "NCBP3", x = colnames(protein.groups.f))
    colnames(protein.groups.f)[NCBP3.cols] <- 
      paste0("Columbia_", colnames(protein.groups.f)[NCBP3.cols])
  }
  
  return(protein.groups.f)
}
