#migration with buffer

#packages
library(stars)
library(tidyr)
library(tidyverse)

#input

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/"

specieslist <- list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species")

rates = read.csv(rates.csv)

results <- lapply(specieslist, function(x) {
  
species = x
rate = rates(x)

#rate1 = rate[1,2]

load(glue::glue("{data}{species}/Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}{species}/Predictions/pred_stars_current.rda"))

#create buffer
rate= 150
buf5 <- rastOut %>% 
  .[.> 0.5] %>% #select pres
  st_as_sf(., coords = c("geometry"), crs = 4326) %>% #as sf
  st_buffer(., dist = rate)  %>% st_union() #buffer with rate

#calculate true pred

true_pred <- pred_stars_list %>%
  .[[1]][,,,1] %>% #select scenario and year
  st_as_sf(., coords = c("geometry"), crs = 4326) %>% #as sf
  true_pred <- st_intersection(., buf) #intersection

})
