require(data.table)

kPathToComplexes <- 'data/corum_complexes_2019-02-22_13_03_24.csv'
kPathToUniprotGeneNames <- "data/uniprot_gene_names.csv"

MakeGroupsList <- function(path.to.complexes) {
  groups.of.complexes <- fread(path.to.complexes)
  output.groups <- list()
  output.subgroup <- list()
  
  for (group in unique(groups.of.complexes$Function)){
    output.groups[[group]] <- 
      groups.of.complexes[groups.of.complexes$Function == group, ]$uniprotID
  }
  
  for (subgroup in unique(groups.of.complexes$ComplexName)){
    output.subgroup[[subgroup]] <- 
      groups.of.complexes[groups.of.complexes$ComplexName == subgroup, ]$uniprotID
  }
  output.groups[['NCBP3']] <- 'Q53F19'
  
  return(list(output.groups, output.subgroup))
}

groups <- MakeGroupsList(kPathToComplexes)
kGroupsList <- groups[[1]]
kSubgroupList <- groups[[2]]
kNoGroupName <- "no group"

kGroupOrder <- fread('data/complex_order.csv')
kGroupOrder <- rbind(kGroupOrder, list('NCBP3', 'NCBP3', 1.1))

kIdColumnName = "Protein IDs"

methods.notnulls <- c('method1','method2','method3','method4','method5','method6',
                     'method7','perseus')

methods.allnulls <- c('method.1','method.2', 'method.3')


kSepMethods <- 
  c("all.mnar", "all.mar", "mnar.0.nonzero", "mnar.1.nonzero",
    "impute.1.col.for.mnar.0.nonzero", "impute.1.col.for.mnar.1.nonzero")

kMethodsMar <- c("mar.stat.1", "mar.stat.2", "mar.stat.3", "mar.stat.4", "mar.stat.5", 
                 "mar.stat.6")
kMethodsMnar <- c('mnar.stat.1', 'mnar.stat.2', 'perseus')

kOldImputationMethods <- 
  list(method1.1=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.1", 
                      method.mnar="mnar.stat.1"),
       method1.2=list(sep.method="impute.1.col.for.mnar.0.nonzero", 
                      method.mar="mar.stat.1", method.mnar="mnar.stat.1"),
       method1.3=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.1", 
                      method.mnar="mnar.stat.2"),
       method2.1=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.2", 
                      method.mnar="mnar.stat.1"),
       method2.2=list(sep.method="impute.1.col.for.mnar.0.nonzero", 
                      method.mar="mar.stat.2", method.mnar="mnar.stat.1"),
       method2.3=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.2", 
                      method.mnar="mnar.stat.2"),
       method3.1=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.3", 
                      method.mnar="mnar.stat.1"),
       method3.2=list(sep.method="impute.1.col.for.mnar.0.nonzero", 
                      method.mar="mar.stat.3", method.mnar="mnar.stat.1"),
       method3.3=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.3", 
                      method.mnar="mnar.stat.2"),
       method4.1=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.4", 
                      method.mnar="mnar.stat.1"),
       method4.2=list(sep.method="impute.1.col.for.mnar.0.nonzero", 
                      method.mar="mar.stat.4", method.mnar="mnar.stat.1"),
       method4.3=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.4", 
                      method.mnar="mnar.stat.2"),
       method5.1=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.5", 
                      method.mnar="mnar.stat.1"),
       method5.2=list(sep.method="impute.1.col.for.mnar.0.nonzero", 
                      method.mar="mar.stat.5", method.mnar="mnar.stat.1"),
       method5.3=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.5", 
                      method.mnar="mnar.stat.2"),
       method6.1=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.6", 
                      method.mnar="mnar.stat.1"),
       method6.2=list(sep.method="impute.1.col.for.mnar.0.nonzero", 
                      method.mar="mar.stat.1", method.mnar="mnar.stat.1"),
       method6.3=list(sep.method="mnar.0.nonzero", method.mar="mar.stat.6", 
                      method.mnar="mnar.stat.2"),
       method7.1=list(sep.method="mnar.1.nonzero", method.mar="mar.stat.6", 
                      method.mnar="mnar.stat.1"),
       method7.2=list(sep.method="impute.1.col.for.mnar.1.nonzero", method.mar="mar.stat.6", 
                        method.mnar="mnar.stat.1"),
       method7.3=list(sep.method="mnar.1.nonzero", method.mar="mar.stat.6", 
                      method.mnar="mnar.stat.2"),
       perseus=list(sep.method="all.mnar", method.mar="", method.mnar="perseus"))

methods.notnulls = unique(str_remove(names(kOldImputationMethods), pattern = "\\..*"))

kColorList <- c(scales::alpha("grey", .3),
               scales::alpha("red", .9),
               scales::alpha("dark green", .9),
               scales::alpha("orange", .9), 
               scales::alpha("turquoise4", .9),
               scales::alpha("green", .9), 
               scales::alpha("blue", .9),
               scales::alpha("black", .9), 
               scales::alpha("brown", .9),
               scales::alpha("violet", .9),
               scales::alpha("cyan", .9))


kColorList2 <- c("red", "blue", "dark green", "orange")

names(kColorList) <- c(kNoGroupName, names(kGroupsList))


kBuffers <- c(18, 12, 7, 20, 10, 14)

kColExp <- c(scales::alpha('red3', 1), scales::alpha('red3', .5), 
             scales::alpha('forestgreen', 1), scales::alpha('forestgreen', .5),
             scales::alpha('royalblue4', 1), scales::alpha('royalblue4', .5),
             scales::alpha('tan1', 1), scales::alpha("tan1", .5))

kAllTargetsPattern <- "Rockefeller|Columbia"
kTargets <- c("Rockefeller_NCBP1", "Rockefeller_NCBP2", "Columbia_NCBP3")

kTargetProteins <- c('Q09161', 'P52298', 'Q53F19')
kGFPprotein <- c('P42212', 'Q9U6Y5')

k.validation.measures <- c('connectivity', 'silhouette', 'dunn', 'BHI_known', 'BHI_GO')
k.validation.measure.borders <- c('min[0,+inf]', 'max[-1,1]', 'max[0,+inf]', 'max[0,1]',
                                  'max[0,1]')
names(k.validation.measure.borders) <- k.validation.measures
k.distance.metrics = c('eLife', 'euclidean')

kN.seeds = 100

