setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")
library(gridExtra)
library(shiny)
library(shinyFiles)

# source("orgaShiny.R")
# ============================================================================
#
#  DESCRIPTION: Data analysis for OrgaMapper
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut f r Molekulare Pharmakologie (FMP)
#               Cellular Imaging - Core facility
#               Campus Berlin-Buch
#               Robert-Roessle-Str. 10
#               13125 Berlin, Germany
#
#         BUGS:
#        NOTES: 
# DEPENDENCIES: ggplot2: install.packages("ggplot2")
#               gridExtra: install.packages("gridExtra")
#
#      VERSION: 0.1.0
#      CREATED: 2020-10-01
#     REVISION: 2020-10-22
#
# ============================================================================
ui <- fluidPage( 
  
  titlePanel("OrgaMapper data analysis"),
  
  sidebarLayout(position = "right",
                
                sidebarPanel(
                  
                  fluidRow(
                    
                    h3("Input:"),
                    
                    shinyDirButton("inputdir", "Input directory", "Upload"),
                    verbatimTextOutput("inputdir", placeholder = TRUE)
                    
                  ),
                  
                  fluidRow(
                    
                    column(width = 5, offset = 0,
                           
                      h3("Processing:"),
                      
                      radioButtons("series_type", 
                                   "Data set type",
                                   choices = list("Singleseries" = 1, 
                                                  "Multiseries" = 2),
                                   selected = 2),
                      
                      textInput(inputId = "series_regex", 
                                label = "Regular expression series number:", 
                                value = "", 
                                width = NULL,
                                placeholder = "(?<=_)\\d*($)"),
                      
                      textInput(inputId = "resultname", 
                                label = "Result Name:", 
                                value = "", 
                                width = NULL,
                                placeholder = "Analysis_test"),

                      sliderInput("feret_limits", 
                                  "Feret's diameter filter:",
                                  min = 0, 
                                  max = 1000, 
                                  value = c(0, 600)),
                      
                      checkboxInput("plot_background_subtract", 
                                    "Subtract background in plots?", 
                                    value = TRUE),
                    ),
                    
                    column(width = 5, offset = 2,
                           
                      h3("Intensity profile:"),
                      
                      checkboxInput("analyze_signal_profiles", 
                                    "Analyze signal profiles?", 
                                    value = FALSE),
                    
                      # plot ranges
                      numericInput("bin_width_norm", 
                                   "Bin width normalized profiles:", 
                                   value = 0.05), 
                      
                      sliderInput("limit_norm", 
                                  "Range of normalized profiles:",
                                  min = 0, 
                                  max = 1, 
                                  value = 0.7),
                      
                      numericInput("bin_width", 
                                   "Bin width raw distance (µm):", 
                                   value = 2),  
                      
                      sliderInput("limit_raw", 
                                  "Range of raw distance plot",
                                  min = 0, 
                                  max = 1000, 
                                  value = 75)
                    )
                  ),
                  
                  fluidRow(
           
                      h3("Processing:"),
                      
                      actionButton("processData", "Plot Data"),
                      
                      h3("Save Results:"),
                      
                      actionButton("saveData", "Save Data")
                      #checkboxInput(inputId = "asCsv", "as .csv", value = FALSE, width = NULL),
                      #checkboxInput(inputId = "asXlsx", "as .xlsx", value = FALSE, width = NULL)
                      
                    ),
                  
                  fluidRow(
                    
                    tags$hr(),
                    tags$p("Author: Christopher Schmied"),
                    tags$p("Contact: schmied@fmp-berlin.de"),
                    tags$p("Cellular Imaging - Core facility")
                    
                  ),
                  
                  fluidRow(
                    tags$img(width = 300, src='Logo.png')
                  )
                  
                ),
                
                mainPanel(
                  
                  tabsetPanel(type = "tabs",
                              # takes outputId plot
                              tabPanel("Cell Measurements", plotOutput("plot", height=1200)),
                              # takes outputId finalTable
                              tabPanel("Table", dataTableOutput(outputId = 'finalTable'))
                              
                  )
                  
                )
                
  )
  
)
server <- function(input, output, session) {
  
  # dir
  shinyDirChoose(
    input, 
    'inputdir', 
    roots = c(home = '~'), 
    filetypes = c('', 'csv'))
  
  global <- reactiveValues(datapath = getwd())
  
  dir <- reactive(input$inputdir)
  
  output$inputdir <- renderText({
    
    global$datapath
    
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$inputdir
               },
               
               handlerExpr = {
                 
                 if ( !"path" %in% names( dir() ) ) 
                   
                   return()
                 home <- normalizePath("~")
                 global$datapath <- file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                 
               })
  
  
  # script
  # gets lists for signal and background
  
  observeEvent(input$processData, {

    inputDirectory <- global$datapath
    
    name_distance = "organelleDistance.csv"
    name_cell_measure = "cellMeasurements.csv"
    
    organelle_distance <- read_collected_files(name_distance, )
    cell_measure <- read_collected_files(name_cell_measure, )

    orga_column = ncol(organelle_distance);
    cell_column = ncol(cell_measure);
    
    
    
    # collects raw data
    table.signal <- collectList(inputDirectory, labelSignal, input$timeRes)
    table.background <- collectList(inputDirectory, labelBackground, input$timeRes)
    
    # calculates average mean intensity per frame
    avg.signal <- calcMean(table.signal)
    avg.background <- calcMean(table.background)
    
    
    
    
    
    # generate final table
    finalTable <- processData(indir, input$frameStim, avg.signal, avg.background)
    
    output$finalTable <- renderDataTable({
      
      finalTable
      
    })
    
    plot.list <- plotMeans(avg.signal, avg.background, finalTable)
    
    
    output$plot <- renderPlot({
      
      plots <- grid.arrange(grobs = plot.list, ncol = 1, top = "Processing results")
      
      print(plots)
      
    })
    
  })
  
  
  observeEvent(input$saveData, {
    
    inputDirectory <- global$datapath
    outputDirectory <- global$datapath
    resultname <- input$resultname
    
    # further settings
    labelSignal = "Spot"
    labelBackground = "background"
    
    inputDirectory <- global$datapath
    
    # collects raw data
    table.signal <- collectList(inputDirectory, labelSignal, input$timeRes)
    table.background <- collectList(inputDirectory, labelBackground, input$timeRes)
    
    # calculates average mean intensity per frame
    avg.signal <- calcMean(table.signal)
    avg.background <- calcMean(table.background)
    
    # generate final table
    finalTable <- processData(indir, input$frameStim, avg.signal, avg.background)
    
    # save files
    
    writeToCsv(outputDirectory, resultname, table.signal, table.background, finalTable)
    
    showNotification("Data saved")
    
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)