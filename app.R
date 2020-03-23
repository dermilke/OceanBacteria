### main script for setting up and running the app ###

### Dependencies ###
library(shiny)
library(shinyWidgets)
library(tidyverse)
library(ggsci)

source("utilities.R")
source("EditMe.R")

### Load Datasets ###
datalist_Atlantic <- list(Count_Data = read_csv(Atlantic_Count_Data),
                          Meta_Data = read_csv(Atlantic_Meta_Data))

datalist_Atlantic$Count_Data <- datalist_Atlantic$Count_Data %>%
  mutate(Richness = 1)

datalist_Pacific <- list(Count_Data = read_csv(Pacific_Count_Data),
                         Meta_Data = read_csv(Pacific_Meta_Data))

datalist_Pacific$Count_Data <- datalist_Pacific$Count_Data %>%
  mutate(Richness = 1)

datalist_Combined <- combineOceanDatalists(datalist_Atlantic, datalist_Pacific)

datalist_Taxonomy <- updateTaxonomy(datalist_Atlantic, datalist_Pacific)

### Setting up Shiny-Server ###

shinyapp(ui, server, options = list(height = 500))