#migration with buffer

#### packages ####
library(stars)
library(tidyr)
library(terra)
library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

#### input ####

#rate = read.csv(rate.csv)

#rate1 = rate[1,2]

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/Chamaenerion angustifolium/"

load(glue::glue("{data}Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}Predictions/pred_stars_current.rda"))

present <- rastOut
future <-  pred_stars_list

#fut126 <- future[[1]] ### 126
#future[[1]][,,,1] ### 126 2026

fut <- future[[1]][,,,1] 

#### Spatial extent ####
wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"

ecoreg <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska


# Transform the CRS of one stars object to match the other
reg_t <- st_transform(reg, st_crs(rastOut))


present <- rastOut %>% st_crop(reg_t)

pres <- rast(present)

#occ
reg_t2 <- st_transform(reg, st_crs(occ_sp_filtered_partition))

occ <- occ_sp_filtered_partition %>% st_crop(reg_t2)

#fut
future <- pred_stars_list[[1]][,,,,1]
reg_t3 <- st_transform(reg, st_crs(future))

fut <- future %>% st_crop(reg_t3)
futsf <- st_as_sf(fut)

#### binary ####

presb <- present[present > 0.5]

#### buffer around presb ####

presb_st <-  st_as_sf(presb, coords = c("geometry"), crs = 4326)
#fut_st <-  st_as_sf(fut, coords = c("geometry"), crs = 4326)

#st_crs(presb_st)$units #is m

buf <- st_buffer(presb_st, dist = 100) 

#st_crs(buf)$units # is 100 m


#### true pred ####

reg_t4 <- st_transform(buf, st_crs(fut))

true_fut <- fut %>% st_crop(reg_t4)

true_fut_0 <- true_fut %>%
  mutate(across(everything(), ~replace(.x, is.na(.x), 0)))

plot(true_fut_0, breaks="equal")
