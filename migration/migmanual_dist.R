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
futsf <- st_as_sf(fut)

#### binary ####

presb <- present[present > 0.5]

presb_st <-  st_as_sf(presb, coords = c("geometry"), crs = 4326)


library(sp)
library(gstat)

xy <- as.data.frame(st_coordinates(presb_st, dims = c("x", "y"))) #extract coords

xyp <- as.data.frame(cbind(xy$X, xy$Y, presb_st$`Chamaenerion angustifolium_MaxEnt_current`))

h <- c("x", "y", "p")
colnames(xyp) = h

xyp$geometry <- st_sfc(lapply(1:nrow(xyp), function(i) st_point(c(xyp$x[i], xyp$y[i]))))

datapoints <- xyp[,-4]

gridded <- expand.grid(x = seq(min(datapoints$x), max(datapoints$x), length.out = 100),
                       y = seq(min(datapoints$y), max(datapoints$y), length.out = 100))

library(phylin)
# Add the interpolated values to the gridded data frame
g <- idw(xyp[-4], gridded)

plot(gridded)


###################

library(gstat)

crs <- 4326
st_as_sf(dat, crs = crs, coords = 
           c("x", "y")) |>
  st_transform(crs) -> dat2


dat = st_as_sf(xyp)

datrc <- st_transform(dat, st_crs(grd))
                         
i <- idw(p~1, dat, grd)

library(stars) |> suppressPackageStartupMessages()
st_bbox(reg_t) |>
  st_as_stars(dx = 10000) |>
  st_crop(reg_t) -> grd



######################
gridded$geometry <- st_sfc(lapply(1:nrow(gridded), function(i) st_point(c(gridded$x[i], gridded$y[i]))))
occs_sf <- st_sf(occs, geometry = ptsb)

j <- left_join(gridded, xyp, by="geometry")

gridded <- as(gridded, "SpatialPixelsDataFrame")
interpolated <- idw(p ~ x, datapoints) 

#### resolution reduction ####

library(sf)
presb_st$geo <- st_round(presb_st$geometry,4)

colnames(presb_st)[1] = "p"
presr <- st_set_precision(presb_st, precision=10^100)


presr <- presr %>%
  group_by(geometry) %>%
  summarize(p = mean(p))

library(dplyr)
pr <- presr %>% distinct(.keep_all = FALSE)

pred <- project(presb_st, res=50000000)

data_1km <- st_transform(presb_st, crs = st_crs(presb_st), 
                         dxy = c(5000000, 5000000))

#### buffer around presb ####

#fut_st <-  st_as_sf(fut, coords = c("geometry"), crs = 4326)

#st_crs(presb_st)$units #is m

#buf <- st_buffer(presb_st, dist = 100) 

#st_crs(buf)$units # is 100 m


#### true pred ####

# reg_t4 <- st_transform(buf, st_crs(fut))
# 
# true_fut <- fut %>% st_crop(reg_t4)
# 
# true_fut_0 <- true_fut %>%
#   mutate(across(everything(), ~replace(.x, is.na(.x), 0)))
# 
# plot(true_fut_0, breaks="equal")

#### distance occs-preds ####

occsf <- st_transform(presb_st, st_crs(futsf))

futsf$dists1 <- rowMins(st_distance(futsf, occsf))

futsf$dists3 <- futsf %>%
            st_distance(., occsf) %>%
            apply(., 1, min, na.rm=TRUE)


# Create a spatial index
occsf_idx <- occsf %>% st_as_sf() %>% st_make_valid() %>% st_intersection()
dist_matrix <- st_distance(futsf, occsf_idx) #nicht alloziieren too big
row_mins <- apply(dist_matrix, 1, min, na.rm = TRUE)

# use lwgeom
library(lwgeom)

dist_matrix <- st_distance(futsf, occsf)
dist_matrix <- st_geod_distance(futsf$geometry, occsf$, tolerance = 0, sparse = FALSE)
futsf$dist <- apply(dist_matrix, 1, min, na.rm = TRUE)

# parallize with future
library(future)
library(future.apply)
plan(multisession)
st_crs(occsf)

occsf <- st_transform(occsf, st_crs(futsf))
future_apply

dist_matrix <- future_apply(futsf, 1, st_distance(futsf, st_transform(occsf, st_crs(futsf))), future.seed = TRUE) #, lazy = FALSE)
row_mins <- apply(dist_matrix, 1, min, na.rm = TRUE)




summary(futsf)
plot(futc_sf$x, futc_sf$dists)

library(ggplot2)
library(viridis)

ggplot(futc_sf, aes(x, y, color = dists)) + 
  geom_point() +
  scale_color_gradientn(colours = viridis::viridis(10)) +
  labs(color = "distance")
