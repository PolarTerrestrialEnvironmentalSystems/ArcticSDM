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

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM"

#region

ecoreg <- st_read(glue::glue("{data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska

#species
specieslist <- as.data.frame(c("Chamaenerion angustifolium", "Achillea millefolium"))
specieslist$rate <- c(50, 150)



yearscen <- as.data.frame(matrix(nrow= 3, ncol=2))
colnames(yearscen) <- c("year", "scen")
yearscen$year <- c(2040, 2070, 2100)
yearscen$scen <- c("spp126", "spp370", "spp858")
  
#for(sp in 1:nrow(specieslist)){

sp = 1
species = specieslist[sp,1]
species = "Achillea millefolium"
load(glue::glue("{data}/Results/{species}/Predictions/pred_stars.rda")) #future
A <- pred_stars[scen][,,,year] 
C <- pred_stars[scen][,,,year]

AC <- c(A[1], C[1], C[1])
AC

library(motif)

#lps
as <- as(AC[1], "SpatRaster")
cs <- as(AC[2], "SpatRaster")
acs <- c(as, cs)
co_occurrence_matrix <- lsp_signature(acs, type = "incoma", neighbourhood = 8)



AC4 <- AC %>%
  mutate(s = .[1] + .[2])


AC[2] + AC[1]
AC4 <- AC[1] %>%
  mutate(s = app(present, present.1, sum))

sum(AC$present.1)


s = rowSums(select(AC, present, present.1))
s                 

AC4 <- AC %>%
  mutate_at(vars(present, present.1), list(s = ~ sum(.)))

AC3 <- AC %>%
  mutate(sum1 = apply(present, present1, sum)

AC3

plot(AC2)
apply(AC[, 1], 2, sum)
AC2
plot(AC2[[1]])
AC2[3]
pred_stars$ssp126
plot(AC2[3])
AC2[[1]][,1]
AC

# Select the two attributes
attr1 <- rast(AC[, 1]) # "present" attribute
attr2 <- rast(AC[, 2])  # "present.1" attribute

# Sum up the attributes for each grid cell
library(terra)
summed_attrs <- attr1 + attr2
extent(attr1

pres <- read_stars(glue::glue("{data}/Results/{species}/{species}_MaxEnt_calibration.tif")) #present
#load(glue::glue("{data}/{species}/tmp/modelTab.rda")) modtab
#load(glue::glue("{data}/{species}/tmp/occurence_thinned.rda"))

 
for(year in 1:2) {
for(scen in 1:2) {
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

#crop future
reg_t3 <- st_transform(reg, st_crs(future))
fut <- future %>% st_crop(reg_t3)
psuit <- as(fut, "SpatRaster")

#### pdist ####

##### prespoints ####
#convert raster to points and only take the presence points

frast <- as(AC[1], "SpatRaster")


occpoints <- terra::as.points(frast, values = TRUE) %>%
          .[which(terra::values(.) > 0.5),]

length(occpoints)
summary(occpoints)

#simeon

binary <- read_stars(glue::glue("{data}/Results/{species}/{species}_MaxEnt_calibration.tif")) %>% 
  setNames("present") %>%
  mutate(present = ifelse(present>0.5, 1, 0)) %>%
  mutate(present = ifelse(is.na(present), 0, present))
A       
C
A <- binary
sum(binary$present)

plot(occpoints)
length(occpoints)
summary(occpoints$`Chamaenerion angustifolium_MaxEnt_calibration.tif`)
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
