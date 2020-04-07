### main script for setting up and running the app ###

### Dependencies ###
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(leaflet)

### Setting up Shiny-Server ###

source("global.R")
source("Source/server.R")
source("Source/ui.R")

shinyApp(ui, server, options = list(height = 500))
