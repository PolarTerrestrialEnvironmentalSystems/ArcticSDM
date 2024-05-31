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

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/Arctic_SDM"

#region

ecoreg <- st_read(glue::glue("{data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska

#speciesb


#get all result species names
#flsl      <- list.dirs(data, 
 #                      recursive = TRUE, full.names = T)

#also unterordner
#need to use structural stuff to extract names
#join with migration_meters.csv

specieslist <- as.data.frame(c("Rubus saxatilis", "Achillea millefolium"))
specieslist$rate <- c(50, 150)


yearscen <- as.data.frame(matrix(nrow= 3, ncol=2))
colnames(yearscen) <- c("year", "scen")
yearscen$year <- c(2040, 2070, 2100)
yearscen$scen <- c("spp126", "spp370", "spp858")
  
#for(sp in 1:nrow(specieslist)){
for(sp in 1:1){

#sp = 1
species = specieslist[sp,1]
load(glue::glue("{data}/Results/{species}/Predictions/pred_stars.rda")) #future
pres <- read_stars(glue::glue("{data}/Results/{species}/{species}_MaxEnt_calibration_restricted_50.tif")) #present
#load(glue::glue("{data}/{species}/tmp/modelTab.rda")) modtab
#load(glue::glue("{data}/{species}/tmp/occurence_thinned.rda"))

 
for(year in 1:3) {

for(scen in 1:3) {
#lapply(1:3, function(scen)
    
#lapply(specieslist, function(scen){
  
scen = 2
year = 3

future <- pred_stars[scen][,,,year] 

#### Spatial extent ####

# crop present
reg_t <- st_transform(reg, st_crs(pres))
present <- pres %>% st_crop(reg_t)
pold <- rast(present)
plot(pold)
#crop future
reg_t3 <- st_transform(reg, st_crs(future))
fut <- future %>% st_crop(reg_t3)
psuit <- as(fut, "SpatRaster")

#### pdist ####

##### prespoints ####
#convert raster to points and only take the presence points

points <- terra::as.points(pold, values = TRUE)
occpoints <- points[which(terra::values(points) > 0.5),]
points
#occpoints <- terra::as.points(pold, values = TRUE) %>%
         # .[which(terra::values(points) > 0.5),]


##### distance ####
#Calculates the distances from each pixel to the nearest presence
CurrDist <- terra::distance(psuit, occpoints)
plot(occpoints)
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

#names
scenn = yearscen[scen, 2]
yearn = yearscen[year, 1]

writeRaster(p, glue::glue("{data}/migResults/{species}_{scenn}_{yearn}.tiff"))


}}}


### Save

library(ggplot2)

dat <- st_as_stars(p)

ggplot() +
  geom_stars(data = dat, show.legend = T)

ggplot() +
  geom_stars(data = dat, show.legend = T) +
  scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
  coord_equal() +
  theme_void()

{
  ggsave(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_predictions.png"), outPredMaps, units = "cm", width = 20, height = 20, bg = "grey70")
  save(pred_stars_list, file = glue::glue("{out_wd}/{sp}/Predictions/pred_stars_list.rda"))  
