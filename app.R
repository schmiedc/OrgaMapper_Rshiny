setwd("/home/schmiedc/FMP_Docs/Repositories/plugins_FMP/orgaMapper_R/")

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
source("plot_intensity_ratio.R")
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
# DEPENDENCIES: shiny: install.packages("shiny")
#               shinyFiles: install.packages("shinyFiles")
#               openxlsx: install.packages("openxlsx")
#               ggplot2: install.packages("ggplot2")
#               gridExtra: install.packages("gridExtra")
#               tidyverse: install.packages("tidyverse")
#               lazyeval: install.packages("lazyeval")
#
#      VERSION: 1.5.0
#      CREATED: 2020-10-01
#     REVISION: 2024-03-28
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
                           
                      h3("General:"),
                      
                      radioButtons("series_type", 
                                   "Data set type",
                                   choices = list("Multi series" = 1,
                                                  "Single series" = 2),
                                   selected = 1),
                      
                      textInput(inputId = "series_regex", 
                                label = "Single series number regex:", 
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

                      checkboxInput("plot_background_subtract", 
                                    "Subtract background in plots?", 
                                    value = TRUE),
                      
                      checkboxInput("analyze_signal_profiles", 
                                    "Map intensity?", 
                                    value = FALSE),
                      
                    ),
                    
                    column(width = 5, offset = 1,
                           
                      h3("Map distance:"),
                      
                      sliderInput("limit_cal_distance", 
                                  "Range of raw distance (µm):",
                                  min = 0, 
                                  max = 1000, 
                                  value = 70),
                      
                      sliderInput("limit_norm_distance", 
                                  "Range of normalized distance:",
                                  min = 0, 
                                  max = 1, 
                                  value = 0.7),
                      
                      h3("Map intensity:"),

                      numericInput("bin_width", 
                                   "Bin width raw distance (µm):", 
                                   value = 2),  
                      
                      sliderInput("limit_raw", 
                                  "Range of raw distance (µm):",
                                  min = 0, 
                                  max = 1000, 
                                  value = 70),
                      
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
      progress = 8
      
      # path to folder where the directories for the measurements are
      directory = "/home/schmiedc/FMP_Docs/Projects/OrgaMapper/2024-02-29_Revision/Tests_manuscript-test/output_bug_test-1/"
      # directory <- global$datapath
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
      
      cal_distance_nucleus <- input$limit_cal_distance
      # norm_distance_nucleus = 0.7
      norm_distance_nucleus <- input$limit_norm_distance
      
      # TODO if file contains series number or the already present column
      # needs to default to something sensible if not possible
      series_type <- input$series_type
      
      single_series = FALSE
      
      if (series_type == 2) {
        
        single_series = TRUE
        
      } else if (series_type == 1) {
        
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
      plots_distance <- file.path(out_dir, "plot_distance_map", fsep = .Platform$file.sep)
      dir.create(plots_distance, showWarnings = FALSE)
      
      # create directory for intensity maps
      plots_intensity <- file.path(out_dir, "plot_intensity_map", fsep = .Platform$file.sep)
      dir.create(plots_intensity, showWarnings = FALSE)
      
      # ==============================================================================
      incProgress(1/progress, detail = paste("Processing distance data. Step: ", 1))
      
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
      
      # Checks if there were measurements in measurement channel
      measureChannelCell = "measureMeanIntensity" %in% colnames(cell_measure)
      measureChannelOrganelle = "measureDetectionPeak" %in% colnames(organelle_distance)
      
      cell_measure_filter <- process_cell_measurements(cell_measure, 
                                                       feret_lower, 
                                                       feret_upper,
                                                       measureChannelCell,
                                                       measureChannelOrganelle,
                                                       feret_filter)
      
      merge_cell_organelle <- process_orga_measurements(cell_measure_filter,
                                                        organelle_distance,
                                                        measureChannelCell,
                                                        measureChannelOrganelle)
      
      merged_summary <- create_summary_table(merge_cell_organelle,
                                             cell_measure_filter)
      
      # ------------------------------------------------------------------------------
      # save processed data
      incProgress(1/progress, detail = paste("Saving distance results. Step: ", 2))
      
      # ------------------------------------------------------------------------------
      # renaming for organelle result tables
      detection_lookup <- c(cell_area = "cellArea",
                       numberOfDetections = "numberDetections",
                       orga_intensity = "orgaMeanIntensity",
                       orga_background = "orgaMeanBackground",
                       x_nucleus_center_mass = "nucleusCenterMassX",
                       y_nucleus_center_mass = "nucleusCenterMassY",
                       measure_intensity = "measureMeanIntensity",
                       measure_background = "measureMeanBackground",
                       orga_intensity_backsub = "orgaMeanIntensityBacksub",
                       measure_intensity_backsub = "measureMeanIntensityBacksub",
                       x_detection = "xDetection",
                       y_detection = "yDetection",
                       orga_distance_pixel = "detectionDistanceRaw",
                       orga_distance_calibrated = "detectionDistanceCalibrated",
                       orga_detection_peak = "orgaDetectionPeak",
                       measure_detection_peak = "measureDetectionPeak",
                       orga_detection_peak_backsub = "orgaDetectionPeakBacksub",
                       measure_detection_peak_backsub = "measureDetectionPeakBacksub",
                       orga_distance_normalized = "detectionDistanceNormalized")
      
        merge_cell_organelle_result <- merge_cell_organelle %>%
          rename(
            any_of(
              
              detection_lookup
              
            )
          )

      # save processed data
      write.xlsx(file = paste0( result_path, "_detection.xlsx", sep = ""), 
                 merge_cell_organelle_result, 
                 sheetName="Sheet1",  
                 colNames=TRUE, 
                 rowNames=TRUE, 
                 append=FALSE, 
                 showNA=TRUE)
      
      # ------------------------------------------------------------------------------
      # renaming for cell results
      cell_lookup <- c(cell_area = "cellArea",
                       orga_numberOfDetections = "numberDetections",
                       orga_intensity = "orgaMeanIntensity",
                       orga_background = "orgaMeanBackground",
                       measure_intensity = "measureMeanIntensity",
                       measure_background = "measureMeanBackground",
                       x_nucleus_center_mass = "nucleusCenterMassX",
                       y_nucleus_center_mass = "nucleusCenterMassY",
                       orga_intensity_backsub = "orgaMeanIntensityBacksub",
                       measure_intensity_backsub = "measureMeanIntensityBacksub",
                       orga_meanDistance_pixel = "detectionDistanceRaw.mean",
                       orga_meanDistance_calibrated = "detectionDistanceCalibrated.mean",
                       measure_intensityOnDetection = "orgaDetectionPeak.mean",
                       orga_intensityOnDetection_backsub = "orgaDetectionPeakBacksub.mean",
                       measure_intensityOnDetection_backsub = "measureDetectionPeakBacksub.mean",
                       orga_meanDistance_normalized = "detectionDistanceNormalized.mean")
      
      merged_summary_result <- merged_summary %>% 
        rename(
          any_of(
            cell_lookup
          )
          )
      
      # merged_summary
      write.xlsx(file = paste0( result_path,  "_cell.xlsx", sep = ""), 
                 merged_summary_result, 
                 sheetName="Sheet1",  
                 colNames=TRUE, 
                 rowNames=TRUE, 
                 append=FALSE, 
                 showNA=TRUE)
      
      # ------------------------------------------------------------------------------
      # plot data
      incProgress(1/progress, detail = paste("Plotting distance maps. Step: ", 3))
      
      cell_plots <- plot_cell_measurements(cell_measure_filter,
                                           plots_distance,
                                           measureChannelCell,
                                           measureChannelOrganelle,
                                           plot_background_subtract)
      
      detection_plots <- plot_detection_measurements(merge_cell_organelle,
                                                     merged_summary,
                                                     plots_distance,
                                                     measureChannelCell,
                                                     measureChannelOrganelle,
                                                     cal_distance_nucleus,
                                                     norm_distance_nucleus,
                                                     plot_background_subtract)
      
      # ------------------------------------------------------------------------------
      incProgress(1/progress, detail = paste("Processing intensity data. Step: ", 4))
      
      output$cell_plots <- renderPlot({
        
        cell <- do.call(grid.arrange, cell_plots)
        
        print(cell)
        
      })
      
      output$orga_plots <- renderPlot({
        
        orga <- do.call(grid.arrange, detection_plots)
        
        print(orga)
        
      })
      
      if (intensity_profiles) {
        
        # ------------------------------------------------------------------------------
        # collect individual files
        print("Computing individual intensity maps")
        individual_intensity_maps <- collect_individual_profiles_new(directory, 
                                                                     series_regex, 
                                                                     single_series, 
                                                                     cell_measure_filter)
        rownames(individual_intensity_maps) <- c()
        head(individual_intensity_maps)
        
        # create intensity ratio data and plots
        
        print("Computing and plotting intensity ratio")
        intensity_ratio_results <- compute_intensity_ration(individual_intensity_maps, 
                                                            10, 
                                                            bin_width, 
                                                            0)
        
        write.xlsx(file = paste0( result_path,  "_intensityRatio.xlsx", sep = ""), 
                   intensity_ratio_results, 
                   sheetName="Sheet1",  
                   colNames=TRUE, 
                   rowNames=TRUE, 
                   append=FALSE, 
                   showNA=TRUE)
        
        plot_intensity_ration(intensity_ratio_results, "orga", plots_intensity)
        
        # group intensity maps
        print("Computing mean of individual intensity maps")
        incProgress(1/progress, detail = paste("Plotting intensity ratio", 5))
        value_lists <- grouped_intensity_map(individual_intensity_maps, plot_background_subtract)
        
        intensity_map_result <- value_lists$raw
        intensity_map_result_norm <- value_lists$norm
        
        head(intensity_map_result)
        # ------------------------------------------------------------------------------
        print("Saving raw intensity maps")
        incProgress(1/progress, detail = paste("Saving intensity data", 6))
        write.xlsx(file = paste0( result_path,  "_intensityProfile.xlsx", sep = ""), 
                   intensity_map_result, 
                   sheetName="Sheet1",  
                   colNames=TRUE, 
                   rowNames=TRUE, 
                   append=FALSE, 
                   showNA=TRUE)
        
        # ------------------------------------------------------------------------------
        print("Plotting intensity maps")
        incProgress(1/progress, detail = paste("Plotting intensity maps", 7))
        
        organelle_profile <- plot_intensity_map(intensity_map_result, 
                                         intensity_map_result_norm, 
                                         "orga", 
                                         bin_width, 
                                         upper_limit,
                                         bin_width_norm,
                                         upper_limit_norm,
                                         plots_intensity)
        
        output$intProfile_organelle <- renderPlot({
          
          intProfile_orga <- do.call(grid.arrange, organelle_profile)
          print(intProfile_orga)
          
        })
        
        # TODO: Fails due to new column length need better way
        if (measureChannelCell && measureChannelOrganelle) {
          
          measure_profile <- plot_intensity_map(intensity_map_result, 
                                              intensity_map_result_norm, 
                                              "measure", 
                                              bin_width, 
                                              upper_limit,
                                              bin_width_norm,
                                              upper_limit_norm,
                                              plots_intensity)
          
          output$intProfile_measure <- renderPlot({
            
            intProfile_meas <- do.call(grid.arrange, measure_profile)
            print(intProfile_meas)
            
          })
          
        }
        
      }
      
      incProgress(1/progress, detail = paste("Done", 8))
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