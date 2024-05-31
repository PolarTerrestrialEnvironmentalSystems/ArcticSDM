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


# Transform the CRS of one stars object to match the other
reg_t <- st_transform(reg, st_crs(rastOut))


present <- rastOut %>% st_crop(reg_t)

pres <- rast(present)

#occ
reg_t2 <- st_transform(reg, st_crs(occ_sp_filtered_partition))

occ <- occ_sp_filtered_partition %>% st_crop(reg_t2)

#fut
future <- pred_stars_list[[1]][,,,1]
reg_t3 <- st_transform(reg, st_crs(future))

fut <- future %>% st_crop(reg_t3)

#### flexsdm ####

futc <- as.data.frame(st_coordinates(fut, dims = c("x", "y"))) #extract coords
futc_xy <- futc[,1:2]

occs <- as.data.frame(st_coordinates(occ, dims = c("x", "y"))) #extract coords
occs$p = 1
colnames(occs)[1] = "x"
colnames(occs)[2] = "y"
#colnames(occs)[3] <- "p"


xp_m <-  extra_eval(
    training_data = occs,#x,y,p
    pr_ab = "p",
    projection_data = futc_xy, #x,y part of real data
    metric = "euclidean",
    univar_comb = TRUE,
    n_cores = 1,
    aggreg_factor = 1)

head(xp_m)

head(futc_sf)


#### truncating #####
#summary(xp_m)
mytrunc_dist <- extra_truncate(
  suit = pres,
  extra = xp_m,
  threshold = 10000,
  trunc_value = 0
)

mytrunc
plot(mytrunc$'10000')
plot(fut)

#### calculate distance manually ####

# predicted points a
ptsa <- st_sfc(lapply(1:nrow(futc), function(i) st_point(c(futc$x[i], futc$y[i]))))
futc_sf <- st_sf(futc, geometry = ptsa)

# occurence points b
ptsb <- st_sfc(lapply(1:nrow(occs), function(i) st_point(c(occs$x[i], occs$y[i]))))
occs_sf <- st_sf(occs, geometry = ptsb)

# Calculate min distances from all points in A to all points in B
library(matrixStats)
futc_sf$dists <- rowMins(st_distance(futc_sf, occs_sf))

xp_m_dist <- xp_m
xp_m_dist$extrapolation <- futc_sf$dists

#### truncating #####
#summary(xp_m)
mytrunc_dist <- extra_truncate(
  suit = pres,
  extra = xp_m_dist,
  threshold = 10000,
  trunc_value = 0
)