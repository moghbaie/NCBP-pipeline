# Mehrnoosh Oghbaie
# 02/10/2019

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
########################################################################################
# Set working directory
########################################################################################

source("Functions.R")
Template <- R6Class("Template",
                    private = list(
                      save_dir = "db/corum_complex_inspected.RData",
                      input_data = "input/max.quant.sign_06282019.rds"
                    ),
                    public = list(
                      complex_enriched = NA,
                      related_complex_enriched = NA,
                      enriched_unique_complexes = NA,
                      average_intensity = NA,
                      related_complex_average_intensity = NA,
                      unrelated_complex_average_intensity = NA,
                      order = data.frame(rbind(
                        c("C complex spliceosome",9),                           
                        c("CTCF-nucleophosmin-PARP-HIS-KPNA-LMNA-TOP complex",14),
                        c("Exon junction complex (EIF4A3, MLN51, MAGOH, Y14)",4),
                        c("Exosome", 10),                                         
                        c("FCP1-associated protein complex", 3),                   
                        c("NEXT complex",11),                                    
                        c("NELF complex (Negative elongation factor complex)",2),
                        c("THO complex",7),                                    
                        c("CBCAP",1),                                           
                        c("Apoptosis and Splicing Associated Protein(ASAP)",5),
                        c("PSAP",6),                                              
                        c("TREX complex",8),                          
                        c("PAXT",12),                                          
                        c("NPC",13)
                      ))
                      
                    ),
                    active = list(
                      complex_protein = function(){
                        load(private$save_dir)
                        return(complex_protein)
                      }
                      ,
                      corum_complex_inspected = function(){
                        load(private$save_dir)
                        return(corum_complex_inspected)
                      }
                      ,
                      max_quant_obj = function(){
                        con <- gzfile("input/max.quant.sign_06282019.rds")
                        data <- readRDS(con)
                        close(con)
                        return(data)
                      },
                      df2 = function(){
                        self$max_quant_obj$GetPvaluesLogFC() %>%
                          dplyr::filter(n.of.anovas.passed>1) %>%
                          dplyr::select(Protein.IDs)
                      }
                      ,
                      df =  function(){
                        dt <- self$max_quant_obj$GetIntensityTable()
                        dt$significant <- ifelse(dt$`Protein IDs` %in% self$df2$Protein.IDs,1,0)
                        dt$uniprotID <- apply(dt, 1, function(x) strsplit(strsplit(x[["Protein IDs"]],";")[[1]][1],"-")[[1]][1])
                        return(dt)
                      },
                      significant_uniprot = function(){
                        return(unlist(unique(self$df$uniprotID[self$df$significant==1])))
                      }, 
                      
                      selected_complex_count_test = function() {
                        load(private$save_dir)
                        return(enrichComplex(self$significant_uniprot,self$complex_protein))
                      }
                      ,
                      anova_result = function(){
                        return(self$max_quant_obj$GetPvaluesLogFC())
                      },
                      conditions = function(){
                        return(self$max_quant_obj$GetMetadataExperiments())
                      }
                    )
)

Template$set("public","enrichProteinPerTarget", function(target){
  ### getting significant proteins in each target protein
  dz <- data.frame(matrix(NA,ncol=6,nrow=0))
  for(tar in unique(self$conditions[[target]])){
    cols <- colnames(self$anova_result)[grepl(tar,colnames(self$anova_result))]
    dz <- rbind(dz,
                self$anova_result %>% dplyr::mutate(pvalue.min = ifelse(do.call(pmin, c(dplyr::select(self$anova_result, matches(paste0("^pvalue.*",tar,".*"))),na.rm=TRUE))<0.05,1,0),
                                                    fold = ifelse(do.call(pmax, c(dplyr::select(self$anova_result,matches(paste0("^log.*",tar,".*"))),na.rm=TRUE))>1,1,0),
                                                    fold.max = do.call(pmax, c(dplyr::select(self$anova_result,matches(paste0("^log.*",tar,".*"))),na.rm=TRUE)),
                                                    target= tar)%>%
                  dplyr::select("Protein.IDs","GENES","fold","pvalue.min","fold.max","target"))
  }
  
  dz$uniprotID <- apply(dz,1, function(x) strsplit(strsplit(x,";")[[1]][1],"-")[[1]][1])
  dz1 <- dz%>% dplyr::filter((fold==1)& (pvalue.min==1))
 
  self$complex_enriched = compareCluster_Complex(dz1,self$selected_complex_count_test,self$complex_protein)
  ##############################################################################################
  ###  filtering list of related complexes
  self$related_complex_enriched = self$selected_complex_count_test%>%
    dplyr::filter(ComplexName %in% self$corum_complex_inspected[["ComplexName"]][self$corum_complex_inspected$Related==1])
  ##############################################################################################
  ### filtering list of related complexes
  self$enriched_unique_complexes = unique(self$complex_enriched[,c("ComplexName","count","complete_members")])
  ##############################################################################################
  ### construct average intensity of proteins per target
  di <- self$max_quant_obj$GetIntensityTable()
  cols2 <- colnames(di)
  for(tar in unique(self$conditions[[target]])){
    di[[paste0("average_",tar)]] <-  round(rowMeans(dplyr::select(di,cols2[grepl(tar,cols2)]),na.rm=TRUE),2)
  }
  di$uniprotID <- apply(di,1, function(x) strsplit(strsplit(x[["Protein IDs"]],";")[[1]][1],"-")[[1]][1])
  
  ### make a table of average intensity in each target per complex
  dt <- data.frame(matrix(NA, nrow=0,ncol=7))
  for(tar in unique(self$conditions[[target]])){
    dx <- data.frame(matrix(NA, ncol=2,nrow=0))
    for(i in 1:dim(self$enriched_unique_complexes)[1]){
      dx  <- rbind(dx ,cbind(rep(as.character(self$enriched_unique_complexes$ComplexName[i]),self$enriched_unique_complexes$count[i]),
                             strsplit(as.character(self$enriched_unique_complexes$complete_members[i]),"\\|")[[1]]))
    }
    dx$target <- tar
    dx$average_intensity <-  di[[paste0("average_",tar)]][match(as.character(dx$V2),di$uniprotID)]
    dx$significant <- ifelse(dz1$target[dz1$target==tar][match(dx$V2,dz1$uniprotID[dz1$target==tar])]==tar,"significant","non-significant")
    res1 <- self$complex_enriched %>% filter(target==tar)
    dx$GeneRatio <- round(res1$GeneRatio[match(as.character(dx$V1),res1$ComplexName)],2)
    dx$p.adj <- round(res1$p.adj[match(as.character(dx$V1),res1$ComplexName)],2)
    dx <-dx[!is.na(dx$V1),]
    dt<- rbind(dt, dx)
  }
  colnames(dt) <- c("ComplexName", "uniprotID", "target", "average_log_intensity", "significant", "GeneRatio","p.adj")
  self$average_intensity = dt
  ##############################################################################################
  ####  average intensity in each target for related complexes
  dt1 <- dt %>% 
    dplyr::mutate(ComplexName = as.character(ComplexName),target = as.character(target)) %>%
    dplyr::group_by(ComplexName,target) %>% 
    dplyr::summarise(average_intensity = mean(average_log_intensity, na.rm=TRUE),GeneRatio=GeneRatio[1],p.adj=p.adj[1])%>%
    dplyr::filter(ComplexName %in% as.character(self$related_complex_enriched$ComplexName), !is.na(p.adj), GeneRatio>=0.5)
  
  dt1 <- data.frame(dt1)
  dt1$order <- as.numeric(as.character(self$order$X2)[match(dt1$ComplexName,as.character(self$order$X1))])
  dt1$ComplexName = factor(dt1$ComplexName, levels=unique(dt1$ComplexName[order(-dt1$order)]), ordered=TRUE)
  self$related_complex_average_intensity = dt1
  #############################################################################################
  #### average intensity in each target for unrelated complexes
  dt2 <- dt %>% 
    dplyr::group_by(ComplexName,target) %>% 
    dplyr::summarise(average_intensity = mean(average_log_intensity, na.rm=TRUE),GeneRatio=GeneRatio[1],p.adj=p.adj[1])%>%
    dplyr::filter(!ComplexName %in% self$related_complex_enriched$ComplexName, !is.na(p.adj), GeneRatio>=0.5)
  self$unrelated_complex_average_intensity = dt2
  
  #plot_related_enriched(self$complex_enriched, self$related_complex_enriched,target) 
  #plot_cluster(dt1,self$related_complex_enriched,title="Related complex enrichment (average intensity& gene ratio)", r=0.0,target)
  #plot_cluster(dt2,self$related_complex_enriched,title="Unrelated complex enrichment(average intensity& gene ratio)", r=0.5,target)
  #plot_cluster(dt2,self$related_complex_enriched,title="Unrelated complex enrichment(average intensity& gene ratio)", r=0.6,target)
})









