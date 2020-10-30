setwd("/data1/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")

packages <- c("shiny", "shinyFiles", "openxlsx", "ggplot2", "gridExtra", "tidyverse", "lazyeval")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

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
#      VERSION: 0.2.0
#      CREATED: 2020-10-01
#     REVISION: 2020-10-28
#
# ============================================================================
ui <- fluidPage( 
  
  titlePanel("OrgaMapper data analysis"),
  
  sidebarLayout(position = "right",
                
                sidebarPanel(
                  
                  fluidRow(
                    
                    h3("Input:"),
                    
                    shinyDirButton("inputdir", "Folder select", "Please select a folder"),
                    verbatimTextOutput("directorypath")
                    
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
                                value = "(?<=_)\\d*($)", 
                                width = NULL,
                                placeholder = "(?<=_)\\d*($)"),
                      
                      textInput(inputId = "resultname", 
                                label = "Result Name:", 
                                value = "Analysis_test", 
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
                                    "Map intensity?", 
                                    value = FALSE),
                 
                      numericInput("bin_width", 
                                   "Bin width raw distance (µm):", 
                                   value = 2),  
                      
                      sliderInput("limit_raw", 
                                  "Range of raw distance (µm):",
                                  min = 0, 
                                  max = 1000, 
                                  value = 75),
                      
                      # plot ranges
                      numericInput("bin_width_norm", 
                                   "Bin width normalized distance:", 
                                   value = 0.05), 
                      
                      sliderInput("limit_norm_profile", 
                                  "Range of normalized distance:",
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
                    
                    column(width = 6, offset = 0,
                      
                      hr(),
                      tags$p("Christopher Schmied"),
                      tags$p("schmied@fmp-berlin.de"),
                      tags$p("Cellular Imaging - Core facility")
                      
                    ),
                    
                    column(width = 4, offset = 0,
                           
                      tags$img(width = 190, src='Logo_mid.png')
                      
                    )
                    
                  )

                ),
                
                mainPanel(  
                  #observeEvent(ignoreNULL = TRUE,
                  
                  tabsetPanel(type = "tabs",
                              # takes outputId plot
                              tabPanel("Cell Measurements", plotOutput("cell_plots", height=1200)),
                              tabPanel("Organelle Measurements", plotOutput("orga_plots", height=1200)),
                              tabPanel("Organelle Intensity Map", plotOutput("intProfile_organelle", height=1200)),
                              tabPanel("Measurement Intensity Map", plotOutput("intProfile_measure", height=1200))
                  )
                  
                )
                
  )
  
)
server <- function(input, output, session) {
  
  volumes <- c(Home = fs::path_home(), 
               "R Installation" = R.home(), 
               getVolumes()())
  
  shinyDirChoose(input, 
                 "inputdir", 
                 roots = volumes, 
                 session = session, 
                 restrictions = system.file(package = "base"))

  global <- reactiveValues(datapath = getwd())

  observe({
    cat("\ninput$directory value:\n\n")
    print(input$inputdir)
  })
  
  output$directorypath <- renderPrint({
    if (is.integer(input$inputdir)) {
      cat("No directory has been selected (shinyDirChoose)")
    } else {
      global$datapath <- parseDirPath(volumes, input$inputdir)
      parseDirPath(volumes, input$inputdir)
    }
  })

  observeEvent(input$processData, {
    
    
    
    tryCatch({
      
      withProgress(message = 'Progress:', value = 0, {
      progress = 7
      
      # path to folder where the directories for the measurements are
      # directory = "/home/schmiedc/Desktop/Test/test_nd2/2020-10-14_output/"
      directory <- global$datapath
      directory <- paste0(directory, .Platform$file.sep)
      
      
      # result_name = "Analysis_test"
      result_name <- input$resultname
      
      # feret_filter = TRUE
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
      intensity_profiles <- input$analyze_signal_profiles
      
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
      incProgress(1/progress, detail = paste("Processing distance data", 1))
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
      incProgress(1/progress, detail = paste("Saving distance results", 2))
      
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
      incProgress(1/progress, detail = paste("Plotting distance maps", 3))
      
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
      
      # ------------------------------------------------------------------------------
      incProgress(1/progress, detail = paste("Processing intensity data", 4))
      
      output$cell_plots <- renderPlot({
        
        cell <- do.call(grid.arrange, cell_plots)
        
        print(cell)
        
      })
      
      output$orga_plots <- renderPlot({
        
        orga <- do.call(grid.arrange, detection_plots)
        
        print(orga)
        
      })

      if (intensity_profiles) {
        
        print("Processing data for intensity maps")
        
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
        print("Saving data of intensity maps")
        value_list <- profile_collected$raw
        value_list_norm <- profile_collected$norm
        
        # delete row names
        rownames(value_list) <- c()
        rownames(value_list_norm) <- c()
        
        # ------------------------------------------------------------------------------
        incProgress(1/progress, detail = paste("Saving intensity data", 5))
        
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
        
        incProgress(1/progress, detail = paste("Plotting intensity maps", 6))
        
        print("Plotting intensity maps")
        organelle_profile <- plot_profiles(value_list, 
                                      value_list_norm, 
                                      "organelle", 
                                      plots_dir, 
                                      plot_background_subtract)

        output$intProfile_organelle <- renderPlot({
          
          intProfile_orga <- do.call(grid.arrange, organelle_profile)
          print(intProfile_orga)
          
        })
        
      
        if (cell_column == 10 && orga_column == 10) {
          
          measure_profile <- plot_profiles(value_list, 
                                            value_list_norm, 
                                            "measure", 
                                            plots_dir, 
                                            plot_background_subtract)
          
          output$intProfile_measure <- renderPlot({
            
            intProfile_meas <- do.call(grid.arrange, measure_profile)
            print(intProfile_meas)
            
          })
           
        } 
  
      } 
      
      incProgress(1/progress, detail = paste("Done", 7))
      }) # end of progess
      
    }, error=function(e) {
      
      message(e)
      showNotification(paste0("WARNING:   ", e), type = 'error')
      
    }, warning=function(w) {
      
      message(w)
      showNotification(paste0("WARNING:   ", w), type = 'warning')
      
    })
      
    })

}

# Run the application
shinyApp(ui = ui, server = server)