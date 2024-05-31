#migmegasdm

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

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/Picea abies/"

load(glue::glue("{data}Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}Predictions/pred_stars_current.rda"))
load(glue::glue("{data}/modelTab_filtered_partition_absence.rda"))
load(glue::glue("{data}/Occurence_thinned.rda"))
present <- rastOut
future <-  pred_stars_list[[2]]

#fut126 <- future[[1]] ### 126
#future[[1]][,,,1] ### 126 2026

fut <- future[[1]][,,,1] 

fut <- future[,,,3] 
fut

#### Spatial extent ####
wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"

ecoreg <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "10583") #choose alaska


# Transform the CRS of one stars object to match the other
reg_t <- st_transform(reg, st_crs(rastOut))


present <- rastOut %>% st_crop(reg_t)
# 
# library(parallel)
# 
# cl <- makeCluster(detectCores())
# presentr <- aggregate(present, fact=100, cores=cl, FUN=mean, by=`Chamaenerion angustifolium_MaxEnt_current`)
# stopCluster(cl)

pres <- rast(present)

#fut
future <- pred_stars_list[[1]][,,,1]
reg_t3 <- st_transform(reg, st_crs(future))

fut <- future %>% st_crop(reg_t3)

futr <- as(fut, "SpatRaster")


######################################################################

#### pdist ####
##### prespoints ####
#convert raster to points and only take the presence points
points <- terra::as.points(pres, values = TRUE)
points2 <- points[which(terra::values(points) > 0.5),]

CurrentPresence <- futr
CurrentPresence[which(terra::values(CurrentPresence) == 0)] <- NA

##### distance ####
#Calculates the distances from each pixel to the nearest presence
CurrDist <- terra::distance(CurrentPresence, points2)
#plot(CurrDist)
#plot(points2, add=T)

##### rate ####

y0 = 1980 #min pres(1980- 2100)
y1 = 2100 #max fut1 
dy = y1-y0
spr = 150 #species specific rate
rate = spr * dy

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
#pdist
#plot(pdist)

#### p ####

pold = pres
psuit = futr

p = pold * pdist * psuit
p1 = pdist * psuit


plot(p1)
plot(pold)
plot(pdist)
plot(psuit)
plot(p)
