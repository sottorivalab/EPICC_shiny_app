library(dplyr)
library(GenomicRanges)
source("functions/atac.R")
source("functions/plotting.R")
source("functions/annotation.R")
source("functions/shinny_app.R")

sessionInfo()

max_extend_by = 10000
n_per_group = 25

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# datasets ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#atac_purity_data = readRDS("data/other/genotyping_estimates_per_sample.rds")
.longer_enh_label = readRDS("data/other/longer_enhancer_labels.rds")

peak_data = readRDS("data/other/peaks.rds")
group_annot = readRDS("data/other/group_annot.rds")
recurrence_summary = readRDS("data/other/recurrence_summary.rds")
insertion_data_annot = readRDS("data/other/insertion_data_annot.rds")

insertion_data = list.files("data/bed_files", "[.]rds$", full.names = TRUE) %>% 
  magrittr::set_names(gsub("[.]rds*", "", basename(.)))

pltA = readRDS("data/other/pltA.rds") #turn_off_clipping()
pltB = readRDS("data/other/pltB.rds") #turn_off_clipping()

# pad labels with spaces do match plot sizes
spaces_to_add = max(nchar(pltB$scales$get_scales("y")$labels)) - max(nchar(pltA$scales$get_scales("y")$labels)) + 10
spaces = paste0(rep(" ", spaces_to_add), collapse="")
scale = pltA$scales$get_scales("y")
pltA = pltA +scale_y_discrete(labels=paste0(spaces, scale$labels), breaks=scale$breaks)

.cna_data = readRDS("data/other/cna_data.rds")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Recurrence heatmap ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

alt_label_tissue = c(
  "Pure" = "Cancer",
  "Adenoma" = "Adenoma",
  "Adenoma (F)" = "Adenoma",
  "Adenoma (G)" = "Adenoma",
  "Adenoma (H)" = "Adenoma"
)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Shinny App ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

slider_row = 
  fluidRow(
    column(3,
           sliderInput(
             "window_width",
             "Context width",
             min = 0,
             max = max_extend_by,
             value = 1000
           )
    ),
    column(3,
           sliderInput(
             "bin_width",
             "Bin width",
             min = 1,
             max = 250,
             value = 100
           )
    )
  )

heatmap_selector_row =
  fluidRow(
    column(8,
      selectInput(
        "heatmap", 
        "Heatmap:",
        c("Drivers" = "drv",
          "Rec. SCCAs" = "rec"),
        )
      ),
    )

heatmap_row =
  fluidRow(
    column(12,
      plotOutput(
        "heatmap",
        click = "plot_click",
        height = "1200px",
        width = "790px"
      )
    )
  )
                                   
track_row = 
  fluidRow(
    column(7,
           plotOutput(
             "plot2",
             height = "300px"
           )
    ),
    column(5,
           plotOutput(
             "plot3",
             height = "300px"
           ),
    )
  )

stats_row = 
  fluidRow(
    column(12,
           verbatimTextOutput("text_stats")
    )
  )

# Define UI for application that draws a histogram
ui =
  fluidPage(
    titlePanel("ATAC-seq Megabulks"),
    heatmap_selector_row,
    heatmap_row,
    slider_row,
    track_row,
    stats_row
  )

# Define server logic required to draw a histogram
server = function(input, output) {
  
  output$heatmap =
    renderPlot(
      if (input$heatmap == "drv") {
        pltA
      } else if (input$heatmap == "rec") {
        pltB
      } else {
        NULL
      }
    )
  
  output$plot2 =
    renderPlot({
      get_click_data(input$plot_click) %>% 
        plot_atac_track(input$window_width, input$bin_width)
    })
  
  output$plot3 =
    renderPlot({
      get_click_data(input$plot_click) %>% 
        plot_sample_stats() 
    })
  
  output$text_stats = 
    renderText({
      get_click_data(input$plot_click) %>% 
        get_stat_text()
    })
  
}
