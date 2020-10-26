setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")
library(gridExtra)
library(shiny)
library(shinyFiles)
library("openxlsx")
library(gridExtra)

source("process_data.R")
source("plot_data.R")
source("process_profiles.R")
source("plot_profiles.R")

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
# DEPENDENCIES: 
#               
#               shiny: install.packages("shiny")
#               shinyFiles: install.packages("shinyFiles")
#               openxlsx: install.packages("openxlsx")
#               ggplot2: install.packages("ggplot2")
#               gridExtra: install.packages("gridExtra")
#               tidyverse: install.packages("tidyverse")
#               lazyeval: install.packages("lazyeval")
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
                           
                      h3("Map distance:"),
                      
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
                      
                      checkboxInput("filter_ferets", 
                                    "Filter by ferets diameter?", 
                                    value = TRUE),
                      
                      sliderInput("limit_norm_distance", 
                                  "Range of normalized distance:",
                                  min = 0, 
                                  max = 1, 
                                  value = 0.7),
                      
                      checkboxInput("plot_background_subtract", 
                                    "Subtract background in plots?", 
                                    value = TRUE),
                    ),
                    
                    column(width = 5, offset = 2,
                           
                      h3("Map intensity:"),
                      
                      checkboxInput("analyze_signal_profiles", 
                                    "Analyze signal profiles?", 
                                    value = FALSE),
                 
                      numericInput("bin_width", 
                                   "Bin width raw distance (µm):", 
                                   value = 2),  
                      
                      sliderInput("limit_raw", 
                                  "Range of raw distance plot",
                                  min = 0, 
                                  max = 1000, 
                                  value = 75),
                      
                      # plot ranges
                      numericInput("bin_width_norm", 
                                   "Bin width normalized profiles:", 
                                   value = 0.05), 
                      
                      sliderInput("limit_norm_profile", 
                                  "Range of normalized profile:",
                                  min = 0, 
                                  max = 1, 
                                  value = 0.7)
                      
                    )
                  ),
                  
                  fluidRow(
           
                      h3("Process & Plot:"),
                      
                      actionButton("processData", "Plot Data"),
                      
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
                              tabPanel("Organelle Measurements", plotOutput("orga_plots", height=1200)),
                              tabPanel("Organelle Profile", plotOutput("intProfile_organelle", height=1200)),
                              tabPanel("Measurement Profile", plotOutput("intProfile_measure", height=1200))
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
    
    # path to folder where the directories for the measurements are
    directory = "/home/schmiedc/Desktop/Test/test_nd2/2020-10-14_output/"
    # directory <- global$datapath
    
    # result_name = "Analysis_test"
    result_name <- input$resultname
    
    # filter for feret's diameter
    feret_filter <- input$filter_ferets
    
    # feret_lower = 0
    feret_lower <- input$feret_limits[1]
    
    # feret_upper = 600
    feret_upper <- input$feret_limits[2]
    
    # determine range for plots
    # norm_distance_nucleus = 0.7
    norm_distance_nucleus <- input$limit_norm_distance
    
    # TODO if file contains series number or the already present column
    # needs to default to something sensible if not possible
    series_type <- input$series_type
    
    single_series = FALSE
    
    if (series_type == 1) {
      
      single_series = TRUE
      
    } else if (series_type == 2) {
      
      single_series = FALSE
      
    }

    # series_regex = "(?<=_)\\d*($)"
    series_regex <- input$series_regex
    
    # plot_background_subtract = TRUE
    plot_background_subtract <- input$plot_background_subtract
    
    # analyze_signal_profiles = TRUE
    analyze_signal_profiles <- input$analyze_signal_profiles
    
    # Binning for intensity profiles
    # upper_limit_norm = 1
    upper_limit_norm <-  input$limit_norm_profile
    # bin_width_norm = 0.05
    bin_width_norm <- input$bin_width_norm
    
    # upper_limit = 75
    upper_limit <- input$limit_raw
    # bin_width = 2
    bin_width <- input$bin_width
    
    # ==============================================================================
    # where to save the data
    out_dir =  directory
    result_path <- file.path(out_dir, result_name, fsep = .Platform$file.sep)
    
    # plot dir
    plots_dir <- file.path(out_dir, "plots", fsep = .Platform$file.sep)
    dir.create(plots_dir, showWarnings = FALSE)
    
    # ==============================================================================
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
    
    cell_column <- ncol(cell_measure)
    orga_column <- ncol(organelle_distance)
    
    cell_measure_filter <- process_cell_measurements(cell_measure, 
                                                     feret_lower, 
                                                     feret_upper,
                                                     cell_column,
                                                     orga_column,
                                                     feret_filter)
    
    merge_cell_organelle <- process_orga_measurements(cell_measure_filter,
                                                      organelle_distance,
                                                      cell_column,
                                                      orga_column)
    
    merged_summary <- create_summary_table(merge_cell_organelle,
                                           cell_measure_filter)
    
    # ------------------------------------------------------------------------------
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
    
    # ------------------------------------------------------------------------------
    # plot data
    cell_plots <- plot_cell_measurements(cell_measure_filter,
                                         plots_dir,
                                         cell_column,
                                         orga_column,
                                         plot_background_subtract)
    
    detection_plots <- plot_detection_measurements(merge_cell_organelle,
                                                   merged_summary,
                                                   plots_dir,
                                                   cell_column,
                                                   orga_column,
                                                   norm_distance_nucleus,
                                                   plot_background_subtract)
    
    if (analyze_signal_profiles) {
      
      name_value_measure = "intDistance.csv"
      
      profile_collected <- collect_individual_profiles(directory,
                                                       name_value_measure,
                                                       cell_measure_filter,
                                                       single_series,
                                                       series_regex,
                                                       upper_limit,
                                                       bin_width,
                                                       upper_limit_norm,
                                                       bin_width_norm)
      
      value_list <- profile_collected$raw
      value_list_norm <- profile_collected$norm
      
      # ------------------------------------------------------------------------------
      write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
                 value_list, 
                 sheetName="Sheet1",  
                 col.names=TRUE, 
                 row.names=TRUE, 
                 append=FALSE, 
                 showNA=TRUE)
      
      write.xlsx(file = paste0( result_path,  "_intensityProfile_norm.xlsx", sep = ""), 
                 value_list_norm, 
                 sheetName="Sheet1",  
                 col.names=TRUE, 
                 row.names=TRUE, 
                 append=FALSE, 
                 showNA=TRUE)
      
      profile_plot <- plot_profiles(value_list, 
                                    value_list_norm, 
                                    "organelle", 
                                    plots_dir, 
                                    plot_background_subtract)
      
      
      if (cell_column == 10 && orga_column == 10) {
        
        measure_profiles <- plot_profiles(value_list, 
                                          value_list_norm, 
                                          "measure", 
                                          plots_dir, 
                                          plot_background_subtract)
        
        
        output$intProfile_measure <- renderPlot({
          
          measure_plots <- do.call(grid.arrange, measure_profiles)
          
          print(measure_plots)
          
        })
        
      } 
      
      output$intProfile_organelle <- renderPlot({
        
        profile_plots <- do.call(grid.arrange, profile_plot)
        
        print(profile_plots)
        
      })

    } 
    
    output$cell_plots <- renderPlot({
      
      cell_plots <- do.call(grid.arrange, cell_plots)
      
      print(cell_plots)
      
    })
    
    output$orga_plots <- renderPlot({
      
      orga_plots <- do.call(grid.arrange, detection_plots)
      
      print(orga_plots)
      
    })
    
  })

  
}

# Run the application
shinyApp(ui = ui, server = server)