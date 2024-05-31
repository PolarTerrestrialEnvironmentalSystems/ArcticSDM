#migloop

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

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species"

#get all result species names
flsl      <- list.dirs(data, 
                       recursive = TRUE, full.names = T)

#also unterordner
#need to use structural stuff to extract names
#join with migration_meters.csv

specieslist <- as.data.frame(c("Chamaenerion angustifolium", "Picea abies"))
specieslist$rate <- c(50, 150)


lapply(specieslist, function(x){

"for" sp in 1:nrow(specieslist)
  
species = specieslist[1,1]
sp = 1

load(glue::glue("{data}/{species}/Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}/{species}/Predictions/pred_stars_current.rda"))
load(glue::glue("{data}/{species}/modelTab_filtered_partition_absence.rda"))
load(glue::glue("{data}/{species}/Occurence_thinned.rda"))

"for" scen in 1:nrow(1:3) year in 1:3
scen = 2
year = 3

future <- pred_stars_list[[scen]][,,,year] #2086

#### Spatial extent ####

#region
wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"

ecoreg <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska

# crop present
reg_t <- st_transform(reg, st_crs(rastOut))
present <- rastOut %>% st_crop(reg_t)
pold <- rast(present)

#crop future
reg_t3 <- st_transform(reg, st_crs(future))
fut <- future %>% st_crop(reg_t3)
psuit <- as(fut, "SpatRaster")

#### pdist ####

##### prespoints ####
#convert raster to points and only take the presence points

points <- terra::as.points(pold, values = TRUE)
occpoints <- points[which(terra::values(points) > 0.5),]

#occpoints <- terra::as.points(pold, values = TRUE) %>%
         # .[which(terra::values(points) > 0.5),]


##### distance ####
#Calculates the distances from each pixel to the nearest presence
CurrDist <- terra::distance(psuit, occpoints)

##### rate ####

#dy = 30
#spr = specieslist[sp,2] #species specific rate
rate = specieslist[sp,2] * 30

##### probability ####
#Creates dispersal probability raster from distance raster

DispersalProbRaster <- function(rate, DistRaster, sh) {
  #Calculates lambda of exponential distribution
  Lambda <- 1 / rate
  
  #When exponential distributions are added, they convolute to a gamma distribution
  GammaProbFun <- function(x) {
    1 - stats::pgamma(x, shape = sh, rate = Lambda)
  }
  
  #Relates distance raster to dispersal probability
  DistProb <- terra::app(DistRaster, fun = function(x){GammaProbFun(x)})
  return(DistProb)
}

pdist <- DispersalProbRaster(rate, CurrDist, sh=1)

#### p ####

p = pdist * psuit

plot(p)



