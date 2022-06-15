library(shiny)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
source("setup.R")

# Run the application 
shinyApp(ui = ui, server = server)

