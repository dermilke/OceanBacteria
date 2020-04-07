source("Source/utilities.R")
source("EditMe.R")

### Load Datasets ###
# Complete Atlantic datalist (Count Data + Meta Data)

datalist_Atlantic <- list(Count_Data = read_csv(Atlantic_Count_Data),
                          Meta_Data = read_csv(Atlantic_Meta_Data)) %>%
  make_depth_groups(.)

datalist_Atlantic$Count_Data <- datalist_Atlantic$Count_Data %>%
  mutate(Richness = 1)

Chloro_datalist_Atlantic <- list(Count_Data = read_csv(Chloro_Atlantic_Count_Data),
                                 Meta_Data = read_csv(Chloro_Atlantic_Meta_Data)) %>%
  make_depth_groups(.)

Chloro_datalist_Atlantic$Count_Data <- Chloro_datalist_Atlantic$Count_Data %>%
  mutate(Richness = 1)

# Complete Pacific datalist (Count Data + Meta Data)
datalist_Pacific <- list(Count_Data = read_csv(Pacific_Count_Data),
                         Meta_Data = read_csv(Pacific_Meta_Data)) %>%
  make_depth_groups(.)

datalist_Pacific$Count_Data <- datalist_Pacific$Count_Data %>%
  mutate(Richness = 1)

Chloro_datalist_Pacific <- list(Count_Data = read_csv(Chloro_Pacific_Count_Data),
                                Meta_Data = read_csv(Chloro_Pacific_Meta_Data)) %>%
  make_depth_groups(.)

Chloro_datalist_Pacific$Count_Data <- Chloro_datalist_Pacific$Count_Data %>%
  mutate(Richness = 1)

# Complete datalist including both Oceans (Count Data + Meta Data)
datalist_Combined <- combineOceanDatalists(datalist_Atlantic, datalist_Pacific)

Chloro_datalist_Combined <- combineOceanDatalists(Chloro_datalist_Atlantic, Chloro_datalist_Pacific)

# Complete Taxonomy of both Oceans (reduced Count Data to enhance computation time)
datalist_Taxonomy <- updateTaxonomy(datalist_Atlantic, datalist_Pacific)

Chloro_datalist_Taxonomy <- updateTaxonomy(Chloro_datalist_Atlantic, Chloro_datalist_Pacific)

# Longhurst-Provinces shapefile
provinces <- rgdal::readOGR("Data/longhurst_v4_2010/Longhurst_world_v4_2010.shp")
province_center <- read_csv("Data/longhurst_v4_2010/ProvinceCenter.csv")
