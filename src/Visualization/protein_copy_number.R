CopyNumber <- function(max.quant.obj){
  intensities <- max.quant.obj$GetIntensityTable()
  case.experiments <- max.quant.obj$GetMetadataExperiments()[max.quant.obj$GetMetadataExperiments()$is.case == T,]$experiment.name
  result.dt <- data.table('experiments' = case.experiments)
  i.s <- c(1, 3, 3)
  j.s <- c(2, 1, 2)
  intensities <- intensities[! intensities$`Protein IDs` == 'P52298-2', ]
  for(k in c(1:length(i.s))){
    i <- i.s[k]
    j <- j.s[k]
    result.dt[[paste0('NCBP', i, '/NCBP', j)]] <- unname(unlist(intensities[grep(kTargetProteins[i], intensities$`Protein IDs`), case.experiments, with = F] / 
                                                                  intensities[grep(kTargetProteins[j], intensities$`Protein IDs`), case.experiments, with = F]))
  }
  metadata <- max.quant.obj$GetMetadataExperiments()
  result.dt <- melt(result.dt, id.vars = "experiments")
  result.dt$location <- str_split_fixed(result.dt$experiments, '_', 2)[, 1]
  result.dt$pretty_name <- metadata[match(result.dt$experiments, metadata$experiment.name), ]$pretty_name
  result.dt$experiments <- paste(result.dt$pretty_name, str_split_fixed(result.dt$experiments, '_', 3)[, 3], sep = '_')
  result.dt[is.na(result.dt$value)]$value <- 0
  p <- ggplot(result.dt, aes(x = experiments, y = log2(value), color = variable)) +
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_color_discrete(name = "") +
    ylab(label = 'log2(iBAQ Intensity ratio)')
  pdf(paste0("out/", MakeDate(), "_copy_number_scatter.pdf"), width = 12, height = 7)
  print(p)
  dev.off()
  
  print("Mann U tets:")
  print(wilcox.test(result.dt[result.dt$variable == "NCBP1/NCBP2"]$value, result.dt[result.dt$variable == "NCBP3/NCBP1" & !is.na(result.dt$value)]$value))
  print(wilcox.test(result.dt[result.dt$variable == "NCBP1/NCBP2"]$value, result.dt[result.dt$variable == "NCBP3/NCBP2" & !is.na(result.dt$value)]$value))
  
  ConfInterval <- function(x){
    error <- qt(0.975, df = length(x) - 1) * sd(x, na.rm = T) / sqrt(length(x))
    return(c(mean(x, na.rm = T) - error, mean(x, na.rm = T) + error))
  }
  means <- list()
  conf.inters <- list()
  for (var in unique(result.dt$variable)){
    means[[var]] <- mean(result.dt[result.dt$variable == var, ]$value, na.rm = T)
    conf.inters[[var]] <- ConfInterval(result.dt[result.dt$variable == var, ]$value)
  }  
  print("Means")
  print(means)
  print("confidence intervals")
  print(conf.inters)
  return(list(means, conf.inters))
}



