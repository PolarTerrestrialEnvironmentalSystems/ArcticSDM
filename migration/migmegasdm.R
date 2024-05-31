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

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/Chamaenerion angustifolium/"

load(glue::glue("{data}Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}Predictions/pred_stars_current.rda"))
load(glue::glue("{data}/modelTab_filtered_partition_absence.rda"))
load(glue::glue("{data}/Occurence_thinned.rda"))
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
# 
# library(parallel)
# 
# cl <- makeCluster(detectCores())
# presentr <- aggregate(present, fact=100, cores=cl, FUN=mean, by=`Chamaenerion angustifolium_MaxEnt_current`)
# stopCluster(cl)

pres <- rast(present)

#occ
reg_t2 <- st_transform(reg, st_crs(occ_sp_filtered_partition))

occ <- occ_sp_filtered_partition %>% st_crop(reg_t2)

occ_sf <- st_as_sf(occ)



#fut
future <- pred_stars_list[[1]][,,,1]
reg_t3 <- st_transform(reg, st_crs(future))

fut <- future %>% st_crop(reg_t3)

futr <- rast(future)
futr <- as(fut, "SpatRaster")
futsf <- st_as_sf(fut)

######################################################################

#convert raster to points and only take the presence points
points <- terra::as.points(pres, values = TRUE)
points2 <- points[which(terra::values(points) > 0.5),]

CurrentPresence <- futr
CurrentPresence[which(terra::values(CurrentPresence) == 0)] <- NA


#Calculates the distances from each pixel to the nearest presence
CurrDist <- terra::distance(CurrentPresence, points2)
plot(CurrDist)
plot(points2, add=T)


dp <- c(CurrDist, CurrentPresence)

dp

names(dp) <- n

n <- c("dist", "p")

names(dp)


#Extends the raster back out to full study area extent
CurrDist <- terra::extend(CurrDist, terra::ext(CurrentPresence))
CurrDist[which(is.na(terra::values(CurrDist)))] <- max(terra::values(CurrDist), na.rm = TRUE) + 1
DistFinal <- terra::mask(CurrDist, CurrentPresence)

#if the data have units that are not degrees or meters, convert distance values to meters.
if(terra::linearUnits(DistFinal) != 0) {
  unitmult <- terra::linearUnits(DistFinal)
  DistFinal <- DistFinal * unitmult
}

#Converts distance (now in meters) to kilometers (for dispersal rate)
DistFinal <- (DistFinal / 1000)
#writes distance raster
# terra::writeRaster(DistFinal,
#             filename = paste0(result_dir, "/", spp, "/", "distance_", Time, "_", Scen, ".grd"),
#             overwrite = TRUE)
rm(CurrentPresence, CurrDist)
gc()
return(DistFinal)
}


plot(DistFinal)

DistRaster = dp

y0 = 1980 #min pres(1980- 2100)
y1 = 2040 #max fut1 (2010 - 2040)
dy = y1-y0
spr = 150 #species specific rate
rate = spr * dy


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
pdist
plot(pdist)

summary(res$lyr.1)
res$lyr.1


p = pold * pdist * psuit

plot(p)
pold = pres
psuit = futr

present <- rastOut %>% st_crop(pdist)
x = seq(1,100,1)
y =1/x

plot(x,y)


#Gets presence pixels in raster 1 but not raster 2
t1nott2 <- function(t1, t2) {
  return(terra::mask(t1, t2, inverse = FALSE, maskvalue = 1, updatevalue = 0))
}

#Calculates overlap between raster 1 and raster 2 presences
overlap <- function(t1, t2) {
  return(terra::mask(t1, t2, inverse = TRUE, maskvalue = 1, updatevalue = 0))
}

#Conducts the actual dispersal rate analyses
FinalDispersal <- function(CurSpp) {
  #Gets species name and relevant dispersal rate
  speciesName = gsub("_", " ", CurSpp)
  if (length(grep(paste0(speciesName), dispersal[, 1])) > 0) {
    #Finds species-specific dispersal rate
    dispRateColumn <- which(unlist(lapply(dispersal, is.numeric)))
    dispersal_rate <- dispersal[grep(paste0(speciesName, "\\s*$"), dispersal[, 1]), dispRateColumn]
    if (length(dispersal_rate) > 1) {
      stop(paste0("More than one dispersal rate found for ", speciesName))
    }
    