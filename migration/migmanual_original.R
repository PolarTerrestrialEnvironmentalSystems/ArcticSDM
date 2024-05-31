#migmanual_original

#### packages ####
library(stars)
library(tidyr)
library(tidyverse)

#### input ####

#rate = read.csv(rate.csv)

#rate1 = rate[1,2]

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/Chamaenerion angustifolium/"

load(glue::glue("{data}Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}Predictions/pred_stars_current.rda"))

fut <- future[[1]][,,,1] 

#input 
present <- rastOut
future <-  pred_stars_list
fut <- future[[1]][,,,1] ### 126 2026

#convert to binary
presb <- present[present > 0.5]

#as sf
presb_st <-  st_as_sf(presb, coords = c("geometry"), crs = 4326)
fut_st <-  st_as_sf(fut, coords = c("geometry"), crs = 4326)

#buffer around presb
buf <- st_buffer(presb_st, dist = 150) 

#select realistic predictions by intersection with the migration radius
true_pred <- st_intersection(fut_st, buf)