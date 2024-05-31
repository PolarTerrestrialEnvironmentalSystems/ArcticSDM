#migcut trial

library(stars)
library(terra)
data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/Results/SDM_Species/Chamaenerion angustifolium/"

load(glue::glue("{data}Predictions/pred_stars_predictions.rda"))
load(glue::glue("{data}Predictions/pred_stars_current.rda"))
load(glue::glue("{data}/Occurence_thinned.rda"))
load(glue::glue("{data}/modelTab_filtered_partition_absence.rda"))

library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"

## Spatial extent
ecoreg <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077")
plot(reg$geometry)

# Transform the CRS of one stars object to match the other
reg_t <- st_transform(reg, st_crs(rastOut))

rastOutsf <- st_as_sf(rastOut)

presentn <- reg_t %>% 
  st_intersection(rastOutsf) %>% st_geometry() %>% st_union() 

fut <- reg_t %>% 
  st_intersection(f) %>% st_geometry() %>% st_union() 

future <- pred_stars_list
occ <- occ_sp_filtered_partition
pres <- rast(present)
#pres
#plot(present)
#plot(occ$geometry, col="red")
# 
# library(ggplot2)
# pa <- ggplot() +
#   geom_stars(data = present) +
#   geom_point(data= modelTab, aes(x=lon, y=lat))
# 
# #pa


#flexsdm
library(flexsdm)
#pr_ab character. Column name with presence and absence (or background points or
                                                               #' pseudo-absences) data (i.e., 1 and 0)
                                                                
#projection_data SpatRaster, data.frame or tibble with environmental condition used for projecting a model (e.g.,
#' a larger, encompassing region, a spatially separate region, or a different time period).
#' If data.frame or tibble is used function will return a tibble object.
#' Otherwise, as SpatRaster object.
#



#pred <- rast(present)
#pred$`Chamaenerion angustifolium_MaxEnt_current`
#Using Mahalanobis distance:


#colnames(modelTab)[1] = "Chamaenerion angustifolium_MaxEnt_current"
modelTab$x <- as.numeric(modelTab$lon)
modelTab$y <- as.numeric(modelTab$lat)
#a =  as.numeric(pred$`Chamaenerion angustifolium_MaxEnt_current`)

occs <- as.data.frame(cbind(modelTab$x, modelTab$y, modelTab$p))
h <- c("x", "y", "p")
colnames(occs) = h
# 
# modelTab <- modelTab[,-1]
# modelTab <- modelTab[,-1]
# modelTab <- modelTab[,-3]
# modelTab <- modelTab[,-5]
# 
# modelTabpred <- modelTab[-3]
# modelTabpred <- modelTab[,1:2]
# 
# 
# modelTabpred = modelTabpred + 10000
# modelTabpred[1:1000,1] = modelTabpred[1:1000,1] + 10000
# #half is far away


pred_stars_list
ssp126 <- pred_stars_list[[1]][,,,2]  #one element from large stars
coords <- st_coordinates(ssp126, dims = c("x", "y")) #extract coords
coords_xy <- coords[,-3]

xp_m <-  extra_eval(
    training_data = occs,#x,y,p
    pr_ab = "p",
    projection_data = coords_xy[1:100,],#x,y part of real data
    metric = "euclidean",
    univar_comb = TRUE,
    n_cores = 1,
    aggreg_factor = 1)


#### manually ####
env_calib2 = occs
env_proj2 = coords_xy[1:2012,]
set <- c(seq(1, nrow(env_proj2), 200), nrow(env_proj2) + 
           1)
rowset <- lapply(seq_len((length(set) - 1)), function(x) {
  set[x]:(set[x + 1] - 1)})
  
v0 <- unique(names(coords_xy))
v0 <- sort(v0)

euc_dist <- function(data, var1, var2){
  if(!is.data.frame(data)) {
    stop('I am so sorry, but this function requires a dataframe\n',
         'You have provided an object of class: ', class(data)[1])
  }
  
  # Make columns of dataframe available w/o quotes
  arguments <- as.list(match.call())
  var1 = eval(arguments$var1, data)
  var2 = eval(arguments$var2, data)
  
  x <- sqrt(sum((var1 - var2)^2))
  return(x)
}

for (i in 1:ncol(env_proj2)) {
envdist <- euc_dist(env_proj2[i, v0], 
                    env_calib2[v0])
envdist <- sapply(data.frame(t(envdist)), min)}

data <- cbind(env_calib2, env_proj2)

colnames(data)[4] = "x2"
envdist <- euc_dist(data, x, x2)
xp_m


env_c = env_calib2[,1:2]

### other functions ####
library(nnspat)
ed <- euc.dist(env_c$x, env_proj2$x)

d <- dist(rbind(env_c, env_proj2))

#### truncating ####
mytrunc <- extra_truncate(
  suit = pres,
  extra = xp_m,
  threshold = 43,
  trunc_value = 0
)

summary(mytrunc)

plot(mytrunc$'43')
#kicks a lot of the data out/sets to 0

#glm_trunc

library(stars)
mytrunc$`0.5`
plot(mytrunc$`0.5`)
plot(present)


######################
#' }
#' 
#extra_eval
training_data = modelTab[1:3,]
pr_ab = "p"
projection_data = modelTabpred
metric = "euclidean" 
univar_comb = FALSE
n_cores = 1
aggreg_factor = 1 

extra_truncate <- function(suit, extra, threshold = 0.3, trunc_value = 0) {
  # names(suit) <- "suit"
  l <- as.list(threshold)
  for (i in 1:length(threshold)) {
    l[[i]] <- suit
    for (ii in 1:terra::nlyr(l[[i]])) {
      l[[i]][[ii]][extra[[1]] > threshold[i]] <- trunc_value
    }
  }
  names(l) <- threshold
  return(l)
}


somevar
pred
10

euc_dist <- function(data, var1, var2){
  if(!is.data.frame(data)) {
    stop('I am so sorry, but this function requires a dataframe\n',
         'You have provided an object of class: ', class(data)[1])
  }
  
  # Make columns of dataframe available w/o quotes
  arguments <- as.list(match.call())
  var1 = eval(arguments$var1, data)
  var2 = eval(arguments$var2, data)
  
  x <- sqrt(sum((var1 - var2)^2))
  return(x)
}
