#### packages ####
library(SSDM)
library(raster)

#### input ####

##### env #####
setwd("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM")
getwd()
files = list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/RCP45/2030",
                   pattern = ".tif$",  
                   full.names = TRUE)

Env0 = stack(files)

#crop alaska from env data

library(rnaturalearth)

usa <- ne_states(country = "United States of America") #select usa
alaska <- subset(usa, usa$iso_3166_2 == "US-AK") #select alaska
Env <- crop(Env0, alaska) #crop env data

##### occurrences #####

Occ0 = read.csv("occurrences/2023-10-23_taxa_climate_df.csv")


Occ <- Occ0 %>%
  filter(Species == unique(Species)[1])

#### SDM ####

SDM_GAM <- modelling('GAM', Occ, 
                     Env, Xcol = 'longitude', Ycol = 'latitude', verbose = FALSE)


#### save model  ####
saveRDS(SDM_GAM, file = "models/SDM_GAM_model.rds")

#### evaluation  ####
knitr::kable(SDM_GAM@evaluation)

plot(SDM_GAM@projection, main = 'SDM\nwith GAM algorithm')


