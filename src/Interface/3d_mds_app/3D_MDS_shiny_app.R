library(plotly)
library(shiny)

source('plot_MDS.R')
max.quant.sign <- readRDS('max_quant_sign.RDS')

list.of.group.names <- unique(max.quant.sign$GetMetadataProteins()$group)
pcoa.df <- CountPcoaDf(max.quant.sign)

ui <- fluidPage(
  titlePanel("Plot settings"),
  sidebarPanel(
    checkboxGroupInput("checkGroups", label = h3("Groups to show"), 
                       choiceNames = list.of.group.names,
                       choiceValues = list.of.group.names,
                       selected = list.of.group.names),
    sliderInput("numberOfAnovas", label = h3("Number of passed ANOVAs:"), min = 1, max = 12,
                value = 1, step = 1, ticks = F),
    width = 3
  ),
  mainPanel(
    plotlyOutput("MDS", height = "800px")
  )
)

server <- function(input, output, session) {
  output$MDS <- renderPlotly({
    Plot3DMDS(max.quant.obj = max.quant.sign, pcoa.df = pcoa.df, 
              groups.to.show = input$checkGroups,
              anovas.passed = input$numberOfAnovas)
  })
}

shinyApp(ui, server)
