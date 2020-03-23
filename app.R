### main script for setting up and running the app ###

### Dependencies ###
library(shiny)
library(shinyWidgets)
library(tidyverse)

source("Source/server.R")
source("Source/ui.R")
source("Source/utilities.R")

source("EditMe.R")

### Load Datasets ###
# Complete Atlantic datalist (Count Data + Meta Data)
datalist_Atlantic <- list(Count_Data = read_csv(Atlantic_Count_Data),
                          Meta_Data = read_csv(Atlantic_Meta_Data))

datalist_Atlantic$Count_Data <- datalist_Atlantic$Count_Data %>%
  mutate(Richness = 1)

# Complete Atlantic datalist (Count Data + Meta Data)
datalist_Pacific <- list(Count_Data = read_csv(Pacific_Count_Data),
                         Meta_Data = read_csv(Pacific_Meta_Data))

datalist_Pacific$Count_Data <- datalist_Pacific$Count_Data %>%
  mutate(Richness = 1)

# Complete datalist including both Oceans (Count Data + Meta Data)
datalist_Combined <- combineOceanDatalists(datalist_Atlantic, datalist_Pacific)

# Complete Taxonomy of both Oceans (reduced Count Data to enhance computation time)
datalist_Taxonomy <- updateTaxonomy(datalist_Atlantic, datalist_Pacific)

### Setting up Shiny-Server ###

shinyapp(ui, server, options = list(height = 500))