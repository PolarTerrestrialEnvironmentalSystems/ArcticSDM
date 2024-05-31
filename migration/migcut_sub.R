#migcut trial

#### packages ####
library(stars)
library(terra)

library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

library(flexsdm)

#### data ####
data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/Chamaenerion angustifolium/"

load(glue::glue("{data}Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}Predictions/pred_stars_current.rda"))
load(glue::glue("{data}/Occurence_thinned.rda"))
load(glue::glue("{data}/modelTab_filtered_partition_absence.rda"))



#### Spatial extent ####
wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"

ecoreg <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska
#plot(reg$geometry)

#presentn
# Transform the CRS of one stars object to match the other
reg_t <- st_transform(reg, st_crs(rastOut))

rastOutsf <- st_as_sf(rastOut)

presentn <- reg_t %>% 
  st_intersection(rastOutsf) %>% st_geometry() %>% st_union() 



# Join with original data to retain attributes
presentnp <- presentn %>%
  left_join(rastoutsf, by = NULL, keep = TRUE) 

presn <- st_as_stars(presentn)
pres <- rast(presn)


#occ
reg_t2 <- st_transform(reg, st_crs(occ_sp_filtered_partition))
occ <- reg_t2 %>% 
  st_intersection(occ_sp_filtered_partition) %>% st_geometry() %>% st_union() 

#fut
future <- pred_stars_list[[1]][,,,1]
futsf <- st_as_sf(future) #still p with 0-1
reg_t3 <- st_transform(reg, st_crs(futsf))
fut <- futsf %>% 
  st_intersection(reg_t3) %>% st_geometry() %>% st_union() #p 0 or 1
#plot(fut)

futsf


#### flexsdm ####

futc <- as.data.frame(st_coordinates(fut, dims = c("x", "y"))) #extract coords
futc_xy <- futc[,1:2]

occs <- as.data.frame(st_coordinates(occ, dims = c("x", "y"))) #extract coords
colnames(occs)[3] <- "p"

xp_m <-  extra_eval(
    training_data = occs,#x,y,p
    pr_ab = "p",
    projection_data = futc_xy, #x,y part of real data
    metric = "euclidean",
    univar_comb = TRUE,
    n_cores = 1,
    aggreg_factor = 1)


#### truncating #####
#summary(xp_m)
mytrunc <- extra_truncate(
  suit = pres,
  extra = xp_m,
  threshold = 100,
  trunc_value = 0
)


# Perform intersection
intersected <- st_intersection(pts, poly) #fut, buf

# Join with original data to retain attributes
result <- mytrunc$'100' %>%
  st_join(futsf, by = NULL, keep = TRUE) 
#Fehler: kann Vektor der Größe 1999.7 GB nicht allozieren


# Join the SpatRaster and data.frame
r_joined <- merge(mytrunc$'100', futsf)

mytrunc$'100'
summary(mytrunc)

plot(mytrunc$'100')

result <- unlist(mytrunc$'100')
result
mytrunc2 <- st_as_sf(result) 
result <- as.data.frame(st_coordinates(mytrunc$'100', dims = c("x", "y"))) #extract coords 

extracted_values <- terra::extract(result, futsf)
#merge mytrunc (0 1) with futsf (0-1)

#mytrunc2 <- mytrunc %>% did not work
  #merge(futsf)


  #mutate(pr = futsf$'2026') did not work

  

#kicks a lot of the data out/sets to 0

#


