setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")
library(gridExtra)
library(shiny)
library(shinyFiles)
library("openxlsx")
library(gridExtra)

source("process_data.R")
source("plot_data.R")

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
                              tabPanel("Cell Measurements", plotOutput("cell_plots", height=1200)),
                              tabPanel("Organelle Measurements", plotOutput("orga_plots", height=1200))
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

    # ==========================================================================
    # Params
    
    # path to folder where the directories for the measurements are
    directory = "/home/schmiedc/Desktop/Test/test_nd2/2020-10-14_output/"
    
    result_name = "Analysis_test"
    
    # filter for feret's diameter
    feret_lower = 0
    feret_upper = 600
    
    # determine range for plots
    norm_distance_nucleus = 0.7
    
    # TODO if file contains series number or the already present column
    # needs to default to something sensible if not possible
    single_series = FALSE
    series_regex = "(?<=_)\\d*($)"
    
    # TODO apply background subtraction for plots
    plot_background_subtract = TRUE
    
    # analyze signal profiles
    analyze_signal_profiles = TRUE
    
    # Binning for intensity profiles
    # or different method for binning
    upper_limit_norm = 1
    bin_width_norm = 0.05
    
    upper_limit = 75
    bin_width = 2
    
    # ==========================================================================
    # where to save the data
    out_dir =  directory
    
    # result_path <- file.path(out_dir, result_name, fsep = .Platform$file.sep)
    
    # plot dir
    
    plots_dir <- file.path(out_dir, plots, fsep = .Platform$file.sep)
    
    dir.create(plots_dir, showWarnings = FALSE)
    
    # ==========================================================================
    name_distance = "organelleDistance.csv"
    name_cell_measure = "cellMeasurements.csv"
    
    organelle_distance <- read_collected_files(directory, 
                                               name_distance, 
                                               single_series, 
                                               series_regex)
    
    
    cell_measure <- read_collected_files(directory, 
                                         name_cell_measure, 
                                         single_series, 
                                         series_regex)
    
    orga_column = ncol(organelle_distance)
    cell_column = ncol(cell_measure)
    
    cell_measure_filter <- process_cell_measurements(cell_measure, 
                                                     feret_lower, 
                                                     feret_upper,
                                                     cell_column,
                                                     orga_column)
    
    merge_cell_organelle <- process_orga_measurements(cell_measure_filter,
                                                      organelle_distance,
                                                      cell_column,
                                                      orga_column)
    
    merged_summary <- create_summary_table(merge_cell_organelle,
                                           cell_measure_filter)
    
    # --------------------------------------------------------------------------
    # save processed data
    write.xlsx(file = paste0( result_path, "_detection.xlsx", sep = ""), 
               merge_cell_organelle, 
               sheetName="Sheet1",  
               col.names=TRUE, 
               row.names=TRUE, 
               append=FALSE, 
               showNA=TRUE)
    
    # merged_summary
    write.xlsx(file = paste0( result_path,  "_cell.xlsx", sep = ""), 
               merged_summary, 
               sheetName="Sheet1",  
               col.names=TRUE, 
               row.names=TRUE, 
               append=FALSE, 
               showNA=TRUE)
    
    # --------------------------------------------------------------------------
    # plot data
    cell_plots <- plot_cell_measurements(cell_measure_filter,
                                         plots_dir,
                                         cell_column,
                                         orga_column)
    
    detection_plots <- plot_detection_measurements(merge_cell_organelle,
                                                   merged_summary,
                                                   plots_dir,
                                                   cell_column,
                                                   orga_column)
    
    
   
    output$cell_plots <- renderPlot({
      
      cell_plots <- do.call(grid.arrange, plot_list_cell)
      
      print(cell_plots)
      
    })
    
    output$orga_plots <- renderPlot({
      
      orga_plots <- do.call(grid.arrange, plot_list_detection)
      
      print(orga_plots)
      
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