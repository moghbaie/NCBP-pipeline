###### Mehrnoosh Oghbaie
### 07/06/2019
### Shiny app

library("data.table")
library("reshape2")
library("dplyr")
library("magrittr")
library("stringr")
library("ggplot2")
library("R6")
library("shiny")
library("plotly")
library("DT")
library("Matrix")
setRepositories(ind = 1,addURLs = c(BioC = "https://bioconductor.org/packages/3.8/bioc"))
library("biomaRt")
source("Functions.R", local=TRUE)
source("Class.R", local=TRUE)

#NCBP <- Template$new()
#NCBP$enrichProteinPerTarget("target")
#save(NCBP, file="db/target.RData")
#NCBP$enrichProteinPerTarget("group")
#save(NCBP, file="db/group.RData")
load(file="db/target.RData")
########################################################
###Make the class


# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("NCBP Significant Protein - complex enrichment"),br(),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select the random distribution type ----
      radioButtons("tar", label = h3("Enrichment level:"),
                   choices = list("Target protein" = "target", "Condition" = "group"), 
                   selected = "target")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                    tabPanel("Related complex bar plot",br(), plotOutput("plot3")),
                  tabPanel("Related complex dot plot",br(), plotlyOutput("plot1")),
                  tabPanel("Unrelated complex dot plot", br(), plotlyOutput("plot2")),
                  tabPanel("Enriched complex table with intensity", 
                           downloadLink("downloadData","Download as csv"),
                           DT::dataTableOutput("mytable"))
      )
      
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression
  d <- reactive({
    #NCBP$enrichProteinPerTarget(input$tar) 
    rm(list=c("NCBP"))
    print(input$tar)
    load(file=paste0("db/",input$tar,".RData"))
    return(NCBP)
  })
 
#dc <- reactive({
#  d()
#})
  
 
  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.
  output$plot1 <- renderPlotly({
    NCBP <- d()
    plotly_cluster(test = NCBP$related_complex_average_intensity,
                   related_complex_enriched = NCBP$related_complex_enriched,
                   title="Related complex enrichment (average intensity& gene ratio)",
                   r=0.0,
                   target=input$tar)
  })
  output$plot2 <- renderPlotly({
    NCBP <- d()
    test <- NCBP$unrelated_complex_average_intensity
    test <- test[test[["GeneRatio"]]>0.6&as.character(test[["ComplexName"]]) %in% as.character(test %>% dplyr::group_by(ComplexName) %>% dplyr::summarize(min(GeneRatio)>0)%>% .$ComplexName) ,] 
    
    test$ComplexName <- as.character(test$ComplexName)
    plotly_cluster(test ,
                   related_complex_enriched = NCBP$related_complex_enriched,
                   title="Unrelated complex enrichment(average intensity& gene ratio)",
                   r=0.6,
                   target=input$tar)
  })
  
  output$plot3 <- renderPlot({
   #target <- input$tar
   NCBP <- d()
    res <- NCBP$complex_enriched
    related_complex_enriched <- NCBP$related_complex_enriched
      result <- data.frame(matrix(NA, nrow=0, ncol=3))
      for(i in unique(res[["target"]])){
        print(i)
        res1 <- res %>% filter(target==i, ComplexName %in% related_complex_enriched$ComplexName)
        x <- cbind(
          as.character(related_complex_enriched$ComplexName),
          res1$GeneRatio[match(related_complex_enriched$ComplexName,res1$ComplexName)],
          rep(i, dim(related_complex_enriched)[1])
        )
        result <- rbind(result,data.frame(x))
      }
      colnames(result) <- c("ComplexName","GeneRatio","target")
      result$GeneRatio <- as.numeric(as.character(result$GeneRatio))
      result$ComplexName <- as.character(result$ComplexName)
      result$GeneRatio[is.na(result$GeneRatio)] <- 0
      ggplot(mapping = aes(x=ComplexName, y=GeneRatio, fill=target))+
        geom_bar(data = related_complex_enriched, fill="gray",stat="identity")+
        geom_bar(data = result, stat="identity",position="dodge")+
        theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14,face="bold"))+
        coord_flip()
 
  }, height = 800)
 
  
  # Generate an HTML table view of the data ----

  output$mytable <- renderDT(

    datatable(d()$average_intensity,
              callback = JS("$('div.dwnld').append($('#downloadData'));"),
              extensions = 'Buttons',
              options = list(
                dom = 'B<"dwnld">frtip',
                buttons = list(
                  "copy"
                )
              )
    )
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data_",gsub(":","-",Sys.time()), ".csv", sep="")
    }, 
    content = function(file) {
      write.csv(data.table(NCBP$average_intensity), file)
    }
  )
  
}
# Create Shiny app ----
shinyApp(ui, server)



