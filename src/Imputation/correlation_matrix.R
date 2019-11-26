require(ggplot2)
require(data.table)
require(ggnewscale)

PlotCorrelationMatrix <- function(intensities, metadata, id.colname = kIdColumnName, 
                                    main = "") {
  intensities <- intensities[, !id.colname, with = F]  
  
  intensities <- as.matrix(intensities)
  intensities <- apply(intensities, 2, as.numeric)
  
  samples.count <- dim(intensities)[[2]]
  
  intensities[is.na(intensities)] <- 0
  
  setorder(metadata, target, is.case, buffer)
  intensities <- intensities[, metadata$experiment.name]
  
  cormat <- round(cor(intensities), 2)
  melted_cormat <- reshape2::melt(cormat)
  
  pdf(paste0("out/", MakeDate(), main, "corr_heatmap.pdf"), paper = "a4")

  case.or.control <- rep(x = "control", length(metadata$is.case))
  case.or.control[metadata$is.case] <- "case"
  
  count.of.replicas <- rle(paste0(metadata$target, "_",
                                  case.or.control, "_", metadata$buffer))
  count.of.exp <- rle(paste0(case.or.control, "_", metadata$target))
  col.exp <- c(scales::alpha('red3',1), scales::alpha('red3',.5), 
               scales::alpha('forestgreen',1), scales::alpha('forestgreen',.5),
               scales::alpha('royalblue4',1), scales::alpha('royalblue4',.5))
  
  exp.end <- cumsum(count.of.exp$lengths) + .5
  exp.start <- c(.5, exp.end[1:(length(exp.end)-1)])
  
  rects <- data.table(xstart = as.numeric(exp.start), 
                      xend = as.numeric(exp.end), 
                      experiment = as.factor(count.of.exp$values))
  
  rects$experiment <- as.factor(rects$experiment)
  
  x.begin = c()
  y.begin = c()
  x.end = c()
  y.end = c()
  
  cur.x = 0.5
  cur.y = 0.5
  
  # making horizontal lines
  for (rep.length in count.of.replicas$lengths) {
    x.begin <- c(x.begin, cur.x)
    y.begin <- c(y.begin, cur.y)
    y.end <- c(y.end, cur.y)
    x.end <- c(x.end, cur.x + rep.length)
    
    cur.y <- cur.y + rep.length
    x.begin <- c(x.begin, cur.x)
    y.begin <- c(y.begin, cur.y)
    cur.x <- cur.x + rep.length
    x.end <- c(x.end, cur.x)
    y.end <- c(y.end, cur.y)
  }
  
  cur.x = 0.5
  cur.y = 0.5
  
  # making vertical lines
  for (rep.length in count.of.replicas$lengths) {
    x.begin <- c(x.begin, cur.x)
    y.begin <- c(y.begin, cur.y)
    y.end <- c(y.end, cur.y + rep.length)
    x.end <- c(x.end, cur.x)
    
    cur.x <- cur.x + rep.length
    x.begin <- c(x.begin, cur.x)
    y.begin <- c(y.begin, cur.y)
    cur.y <- cur.y + rep.length
    x.end <- c(x.end, cur.x)
    y.end <- c(y.end, cur.y)
  }
  
  my.lines <- data.frame(x = x.begin, y = y.begin, xend = x.end, yend = y.end)
  
  grid.h <-
    data.frame(
      x = c(rep(0.5, samples.count)),
      y = c(1:samples.count) + 0.5,
      xend = c(rep(samples.count + 0.5, samples.count)),
      yend = c(1:samples.count) + 0.5
    )
  grid.v <-
    data.frame(
      y = c(rep(0.5, samples.count)),
      x = c(1:samples.count) + 0.5,
      yend = c(rep(samples.count + 0.5, samples.count)),
      xend = c(1:samples.count) + 0.5
    )
  
  g <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white", size = 0.3) +
    scale_fill_gradient2(
      low = "skyblue3",
      high = "firebrick2",
      mid = "white",
      midpoint = 0.5,
      limit = c(0, 1),
      space = "Lab",
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    coord_fixed() +
    geom_segment(
      data = grid.h,
      aes(x, y, xend = xend, yend = yend),
      size = 0.2,
      inherit.aes = F,
      color = "palegreen4"
    ) +
    geom_segment(
      data = grid.v,
      aes(x, y, xend = xend, yend = yend),
      size = 0.2,
      inherit.aes = F,
      color = "palegreen4"
    ) +
    geom_segment(
      data = my.lines,
      aes(x, y, xend = xend, yend = yend),
      size = 0.5,
      inherit.aes = F, color = "black"
    ) + labs(x="", y="")
  
  
  g <- g +
    geom_rect(data = rects, mapping = aes(xmin = xstart, xmax = xend, 
                                          ymin = -3, ymax = 0),
              fill = col.exp, inherit.aes = F) +
    geom_rect(data = rects, mapping = aes(xmin = -3, xmax = 0,
                                          ymin = xstart, ymax = xend),
              fill = col.exp, inherit.aes = F, show.legend = T)
  
  print(g)
  plot.new()
  legend("topright", legend=rects$experiment, col = col.exp, lty = 1, lwd=5)
  
  
  dev.off()
  
}

CheckExperimentsCorr <- function(path = "./data/txt_noGA/proteinGroups.txt", pattern,
                                 sep.method, method.mar, method.mnar) {
  max.quant <- BuildMaxQuant(path, pattern = pattern,
                             groups.list = kGroupsList, subgroup.list = kSubgroupList)
  max.quant$MakePrettyExperimentNames()
  max.quant$
    RemoveContaminantsReversed()$
    SumMAGOH(save.result.table = F)$
    LogTricky() 
  
  PlotCorrelationMatrix(intensities = max.quant$GetIntensityTable(),
                        metadata = max.quant$GetMetadataExperiments(), 
                        main = "_LFQ_before_imputation")
  
  max.quant$ImputeAll(id.column = kIdColumnName, sep.method = sep.method, 
                      method.mar = method.mar, method.mnar = method.mnar, sigma.const = 4,
                      remove.outliers = T,
                      save.result.table = F, seed = 46)
  
  PlotCorrelationMatrix(intensities = max.quant$GetIntensityTable(),
                        metadata = max.quant$GetMetadataExperiments(), main = "_LFQ_after_imputation")
  
}