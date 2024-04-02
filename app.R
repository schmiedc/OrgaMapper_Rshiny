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
      # directory = "/home/schmiedc/FMP_Docs/Projects/OrgaMapper/2024-02-29_Revision/Tests_manuscript-test/output_bug_test-1/"
      # directory = "/home/schmiedc/FMP_Docs/Projects/OrgaMapper/2024-02-29_Revision/Feature_tests/output_single_measure_membrane_1/"
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
      
      # ========================================================================
      # where to save the data
      out_dir =  directory
      result_path <- file.path(out_dir, result_name, fsep = .Platform$file.sep)
      
      # plot dir
      plots_distance <- file.path(out_dir, "plot_distance_map", fsep = .Platform$file.sep)
      dir.create(plots_distance, showWarnings = FALSE)
      
      # create directory for intensity maps
      plots_intensity <- file.path(out_dir, "plot_intensity_map", fsep = .Platform$file.sep)
      dir.create(plots_intensity, showWarnings = FALSE)
      
      # ========================================================================
      # Data processing 
      # ========================================================================
      incProgress(1/progress, detail = paste("Processing distance data. Step: ", 1))
      
      name_distance = "organelleDistance.csv"
      name_cell_measure = "cellMeasurements.csv"
      
      cell_measure <- read_collected_files(directory, 
                                           name_cell_measure, 
                                           single_series, 
                                           series_regex)
      
      organelle_distance <- read_collected_files(directory, 
                                                 name_distance, 
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
      
      
      # ------------------------------------------------------------------------
      # check if measurements from membrane exist
      name_distance_membrane = "organelleDistanceFromMembrane.csv"
      
      distance_membrane_file_exists = file.exists(paste0(directory, name_distance_membrane))
      
      if (distance_membrane_file_exists) {
        
        cat(file=stderr(), "Distance from membrane exists")

        # get files for organelle distance from membrane
        organelle_distance_membrane <- read_collected_files(directory,
                                                            name_distance_membrane,
                                                            single_series,
                                                            series_regex)
        
        merge_cell_organelle_membrane <- process_orga_measurements(cell_measure_filter,
                                                          organelle_distance_membrane,
                                                          FALSE,
                                                          FALSE)
        
        merged_summary_membrane <- create_summary_table(merge_cell_organelle_membrane,
                                               cell_measure_filter)
        
        # filter out redundant columns
        keep_columns_detection_membrane <- c("identifier", 
                          "series", 
                          "cell", 
                          "detection", 
                          "detectionDistanceRaw", 
                          "detectionDistanceCalibrated",
                          "detectionDistanceNormalized")
        
        merge_cell_organelle_membrane_2 <- merge_cell_organelle_membrane[, keep_columns_detection_membrane, drop = FALSE]
        
        # rename the column headers for membrane distance
        merge_cell_organelle_membrane_2 <- merge_cell_organelle_membrane_2 %>% rename(
          orga_meanDistance_membrane_pixel = detectionDistanceRaw,
          orga_meanDistance_membrane_calibrated = detectionDistanceCalibrated,
          orga_meanDistance_membrane_normalized = detectionDistanceNormalized)
        
        # filter out redundant columns
        keep_columns_detection_summary <- c("identifier", 
                                             "series", 
                                             "cell",
                                             "detectionDistanceRaw.mean", 
                                             "detectionDistanceCalibrated.mean",
                                             "detectionDistanceNormalized.mean")
        
        merged_summary_membrane_2 <- merged_summary_membrane[, keep_columns_detection_summary, drop = FALSE]
        
        # rename the column headers for membrane distance
        merged_summary_membrane_2 <- merged_summary_membrane_2 %>% rename(
          orga_meanDistance_membrane_pixel = detectionDistanceRaw.mean,
          orga_meanDistance_membrane_calibrated = detectionDistanceCalibrated.mean,
          orga_meanDistance_membrane_normalized = detectionDistanceNormalized.mean)
        
        
        # merge into table that contain all detections
        merge_cell_organelle <- merge(merge_cell_organelle,
                                      merge_cell_organelle_membrane_2,
                                      by = c("identifier", "series", "cell", "detection"))
        
        # merge into table that contain summaries
        merged_summary <- merge(merged_summary,
                                merged_summary_membrane_2,
                                by = c("identifier", "series", "cell"))
        
      } else {
        
        cat(file=stderr(), "Distance from membrane does not exist")
     
      }
      
      # ========================================================================
      # save processed data
      # ========================================================================
      incProgress(1/progress, detail = paste("Saving distance results. Step: ", 2))
      
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
                       orga_distance_nucleus_pixel = "detectionDistanceRaw",
                       orga_distance_nucleus_calibrated = "detectionDistanceCalibrated",
                       orga_detection_peak = "orgaDetectionPeak",
                       measure_detection_peak = "measureDetectionPeak",
                       orga_detection_peak_backsub = "orgaDetectionPeakBacksub",
                       measure_detection_peak_backsub = "measureDetectionPeakBacksub",
                       orga_distance_nucleus_normalized = "detectionDistanceNormalized")
      
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
                       orga_meanDistance_nucleus_pixel = "detectionDistanceRaw.mean",
                       orga_meanDistance_nucleus_calibrated = "detectionDistanceCalibrated.mean",
                       measure_intensityOnDetection = "orgaDetectionPeak.mean",
                       orga_intensityOnDetection_backsub = "orgaDetectionPeakBacksub.mean",
                       measure_intensityOnDetection_backsub = "measureDetectionPeakBacksub.mean",
                       orga_meanDistance_nucleus_normalized = "detectionDistanceNormalized.mean")
      
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
      
      # ========================================================================
      # plot distance data
      # ========================================================================
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
      
      # ------------------------------------------------------------------------
      incProgress(1/progress, detail = paste("Display plots. Step: ", 4))
      
      output$cell_plots <- renderPlot({
        
        cell <- do.call(grid.arrange, cell_plots)
        
        print(cell)
        
      })
      
      output$orga_plots <- renderPlot({
        
        orga <- do.call(grid.arrange, detection_plots)
        
        print(orga)
        
      })

      # ========================================================================
      # Process intensity profiles 
      # ========================================================================
      
      if (intensity_profiles) {
        incProgress(1/progress, detail = paste("Processing intensity profiles", 5))
        
        intensityProfile_nucleus = "intensityDistance.csv"

        # collect individual files
        print("Collecting individual intensity maps")
        individual_intensity_maps <- collect_individual_profiles_new(directory, 
                                                                     series_regex, 
                                                                     single_series,
                                                                     intensityProfile_nucleus,
                                                                     cell_measure_filter,
                                                                     measureChannelCell)
        rownames(individual_intensity_maps) <- c()
        
        # Identify if a measurement channel is present
        measureChannelIntensity = "mean_measureIntensity" %in% colnames(individual_intensity_maps)
        
        # ----------------------------------------------------------------------
        # create intensity ratio data and plots
        print("Computing intensity ratio")
        intensity_ratio_results <- compute_intensity_ratio(individual_intensity_maps, 
                                                           10, 
                                                           bin_width, 
                                                           0)
        
        # ----------------------------------------------------------------------
        print("Computing mean of individual intensity maps")

        value_lists <- grouped_intensity_map(individual_intensity_maps, 
                                             plot_background_subtract,
                                             measureChannelCell,
                                             measureChannelIntensity)
        
        intensity_map_result <- value_lists$raw
        intensity_map_result_norm <- value_lists$norm
        
        # ----------------------------------------------------------------------
        # TODO: Compute mean intensity for membrane measurements
        # check if measurements from membrane exist
        name_distance_membrane = "organelleDistanceFromMembrane.csv"
        
        distance_membrane_file_exists = file.exists(paste0(directory, name_distance_membrane))
        # distance_membrane_file_exists = FALSE 
        
        if (distance_membrane_file_exists) {
          
          cat(file=stderr(), "Distance from membrane exists")
          
          intensityProfile_membrane = "intensityDistanceFromMembrane.csv"
          
          
          
          # collect individual files
          print("Collecting individual intensity maps for membrane measurement")
          individual_intensity_maps_membrane <- collect_individual_profiles_new(directory,
                                                                                series_regex,
                                                                                single_series,
                                                                                intensityProfile_membrane,
                                                                                cell_measure_filter,
                                                                                measureChannelCell)
          
          measureChannelIntensity_membrane = "mean_measureIntensity" %in% colnames(individual_intensity_maps_membrane)
          
          grouped_intensity_maps_membrane <- grouped_intensity_map(individual_intensity_maps_membrane, 
                                                                   plot_background_subtract,
                                                                   measureChannelCell,
                                                                   measureChannelIntensity_membrane)
          
          intensity_map_membrane_result <- grouped_intensity_maps_membrane$raw
          intensity_map_membrane_result_norm <-  grouped_intensity_maps_membrane$norm
          
          write.xlsx(file = paste0( result_path,  "_intensityProfile_Membrane.xlsx", sep = ""), 
                     intensity_map_membrane_result, 
                     sheetName="Sheet1",  
                     colNames=TRUE, 
                     rowNames=TRUE, 
                     append=FALSE, 
                     showNA=TRUE)
          
        } else {
          
          cat(file=stderr(), "Distance from membrane does not exist")
          
        }
        # ----------------------------------------------------------------------
        print("Plotting intensity data")
        plot_intensity_ratio(intensity_ratio_results, "orga", plots_intensity)
        
        incProgress(1/progress, detail = paste("Plotting intensity maps", 6))
        
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
        
        # Intensity maps based on measurement channel
        if (measureChannelCell && measureChannelIntensity) {
          
          print("Plotting intensity maps for measure channel")
          
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
        
        # ----------------------------------------------------------------------
        print("Saving intensity data")
        incProgress(1/progress, detail = paste("Saving intensity data", 7))
        
        write.xlsx(file = paste0( result_path,  "_intensityProfile_Nucleus.xlsx", sep = ""), 
                   intensity_map_result, 
                   sheetName="Sheet1",  
                   colNames=TRUE, 
                   rowNames=TRUE, 
                   append=FALSE, 
                   showNA=TRUE)
        
        write.xlsx(file = paste0( result_path,  "_intensityRatio_Nucleus.xlsx", sep = ""), 
                   intensity_ratio_results, 
                   sheetName="Sheet1",  
                   colNames=TRUE, 
                   rowNames=TRUE, 
                   append=FALSE, 
                   showNA=TRUE)
        
      }
      # ========================================================================
      # End processing
      # ========================================================================
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