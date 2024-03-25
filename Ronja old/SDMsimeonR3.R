
#### Preparation  ####

library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)


 wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
#  wd <- "/Volumes/projects/bioing/data/ArcticSDM/"
# wd <- '/bioing/data/ArcticSDM/'
# map_data <- "~/Documents/ArcticSDM_data/"
# map_data <- '/bioing/data/ArcticSDM/data/'
 map_data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/data/"
 
##### Projection ####
proj <- "+proj=laea +lon_0=-170 +lat_0=90"
  
##### Maps ####
ecoreg   <- st_read(glue::glue("{map_data}Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

library(rnaturalearth)
 
map_wrld <- ne_countries(scale = 50, returnclass = 'sf')  %>% 
  st_intersection(ecoreg %>% st_union()) %>% dplyr::select(c('name','admin'))

  #st_read(glue::glue("{map_data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
#path did not work
#need to download data
#alaska is not possible with the downloaded data
#would need to use other approach
#usa <- ne_states(country = "United States of America") #select usa
#alaska <- subset(usa, usa$iso_3166_2 == "US-AK") #select alaska

##### Subset ####
map <- map_wrld %>% filter(name=="Canada")  %>% 
  st_transform(proj)

##### Resolution #####
res <- 5000
grd_map <- st_rasterize(map, st_as_stars(st_bbox(map), 
             dx = res, dy = res, 
             values = NA_real_)) %>% setNames("grd_map")


##################
#### GBIF Occ ####
##################

gbif_grid <- st_read(glue::glue("{map_data}occurance_grid.shp"))

ids <- which(c(gbif_grid %>% st_transform(st_crs(map)) %>%
  st_intersects(map, sparse = FALSE)))

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = gbif_grid[ids,] %>% st_transform(st_crs(map)), fill = NA)
 
##### Read in data files #####
gbif_list <- tibble(fls = list.files(glue::glue("{wd}/data/GBIF/Plantea"))) %>%
  mutate(ids = as.numeric(gsub("_gbif_plantea.csv", "", fls)))

##### Filter #####
#mcapply is a parallelized version of lapply, 
#it returns a list of the same length as X, 
#each element of which is the result of applying FUN to the corresponding element of X. 
#It relies on forking and hence is not available on Windows unless mc.cores = 1

#cannot select kingdom etc for gbif_list[1]
#choose only [2]
#quicker calculating for less data

gbifTab0    <- parallel::mclapply(which(gbif_list$ids%in%ids), function(x) {
  read_csv(glue::glue("{map_data}/GBIF/Plantea/{gbif_list$fls[x]}"), progress = F, show_col_types = FALSE) %>%
    dplyr::select(kingdom, phylum, class, order, family, genus, species, scientificName, countryCode, 
                  decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters) %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species)) 
  }, mc.cores = 1) %>% Reduce("rbind",.) %>%
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
  st_transform(st_crs(map))

#with lapply and only csv 2
library(dplyr)
gbifTab    <- lapply(which(gbif_list$ids%in%ids[2:2]), function(x) {
  read_csv(glue::glue("{map_data}/GBIF/Plantea/{gbif_list$fls[x]}"), progress = F, show_col_types = FALSE)  %>%
  dplyr::select(kingdom, phylum, class, order, family, genus, species, scientificName, countryCode, 
                decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species)) 
}) %>% Reduce("rbind",.) %>%
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
  st_transform(st_crs(map))


##### Thinning #####
res_thin <- 15000
thin_map <- st_rasterize(map, st_as_stars(st_bbox(map), 
                                         dx = res_thin, dy = res_thin, 
                                         values = NA_real_)) %>% setNames("grd_map")

#select species with 10at least 50 records
gbif_thin <- (gbifTab %>% group_split(species))[sapply(gbifTab %>% group_split(species), nrow)>50] %>%
  parallel::mclapply(function(x) {
    st_rasterize(x %>% dplyr::select(species), thin_map) %>% st_as_sf(na.rm = FALSE) %>%
      mutate(poly = st_as_sf(thin_map, na.rm = FALSE) %>% pull(grd_map)) %>%
      st_centroid() %>%
      rownames_to_column(var = "id") %>%
      filter(!is.na(poly) & !is.na(ID) & ID > 1) %>%
      mutate(species = x$species[1]) %>%
      dplyr::select(id, species) %>% suppressWarnings()
    }, mc.cores = 1) %>% Reduce("rbind", .)

write_csv(gbif_thin, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/data/GBIF/gbif_thin_test2.csv")

############################
#### Environmental data ####
############################

#gbif_thin <- read_csv("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/data/GBIF/gbif_thin_test2.csv")
tt <- gbif_thin %>% dplyr::select(id) %>% filter(!duplicated(id))

env_list  <- list.files(glue::glue("{wd}environment/calibration"), pattern = ".tiff")

env_stars <- lapply(env_list, function(x) {
  st <- read_stars(glue::glue("{wd}environment/calibration/{x}"))
  st_extract(st, tt %>% st_transform(st_crs(st))) %>% st_as_sf() %>% st_drop_geometry()
}) %>% Reduce("cbind", .)
#works with produced gbif_thin but not if read in 
#only select env data with in the grid cells we want (like occ, matching extent)

##################
#### Species Loop ####
##################

sp <- gbif_thin %>% st_drop_geometry() %>% group_by(species) %>% summarise(sample = n()) %>%
  arrange(desc(sample))

s <- sp %>% slice(1) %>% pull(species) #select first (most abundant) species

#generating absence data
#more columns, less rows, so one column for each species  and rows for grid cell
#present species get 1, 
#not mentioned species at position get NA
#sum up, order
#nr of abs should not be longer than nr of pres

occ_wide   <- gbif_thin %>% st_drop_geometry() %>% mutate(vals = 1) %>% pivot_wider(id_cols = id, names_from = species, values_from = vals)
id_occ     <- occ_wide[which(!is.na(occ_wide %>% pull(s))),] %>% pull(id)
id_abs_tab <- occ_wide[-which(!is.na(occ_wide %>% pull(s))),] %>% 
  mutate(sum = apply(occ_wide[-which(!is.na(occ_wide %>% pull(s))),-1], 1, sum, na.rm = T)) %>%
  arrange(desc(sum))
id_abs     <- id_abs_tab[{if(nrow(id_abs_tab)<length(id_occ)) 1:nrow(id_abs_tab) else 1:length(id_occ)},] %>% pull(id)

mod_tab <- gbif_thin[gbif_thin$id%in%id_occ,] %>% mutate(p = 1) %>%
  bind_rows(gbif_thin[gbif_thin$id%in%id_abs,] %>% mutate(p = 0)) %>% st_drop_geometry() %>%
  dplyr::select(id, p) %>% right_join(env_stars %>% mutate(id=tt$id), by = 'id') %>%
  dplyr::select(-id)
  
write_csv(mod_tab, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/data/GBIF/mod_tab_2_Rubus_odoratus.csv")

#### SDM ####

mod_tab <- read_csv("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/data/GBIF/mod_tab_2_Rubus_odoratus.csv")

##### GLM #####
formula <- paste0("p ~ ", paste(colnames(mod_tab)[-1], collapse = " + "))
glm_mod <- glm(formula, family = binomial(link = "logit"), data = mod_tab)


##### GAM #####
library(mgcv)
#gam_mod2 <- gam(formula, family = binomial(link = "logit"), data = mod_tab)
formula
gam_mod <- gam(p ~ p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14 + p14.1 + p16 + p17 + p18 + p19 + wc2.1_2.5m_prec_2015.tiff + wc2.1_2.5m_tmax_2015.tiff + wc2.1_2.5m_tmin_2015.tiff,
    family = binomial(link = "logit"), data=mod_tab)  

##### MARS #####
library(earth)
mars_mod <- earth(p ~ .,  data = mod_tab)  

##### MAXENT #####

library(dismo)
# set up the argument list
args_list <- c('responsecurves=TRUE',
               'jackknife=TRUE',
               'pictures=TRUE',
               'autofeature=FALSE',
               'linear=TRUE',
               'quadratic=TRUE',
               'product=TRUE',
               'threshold=TRUE',
               'hinge=TRUE')
               #paste0('betamultiplier=', rm))

maxent_mod <- maxent(x=mod_tab[-1] , p=mod_tab$p, a=NULL, removeDuplicates=TRUE, nbg=10000,
             path=wd, args=args_list)


#### Predictions ####
res_pred <- 2500
grd_out  <- st_rasterize(map, st_as_stars(st_bbox(map), 
                                         dx = res_pred, dy = res_pred, 
                                         values = NA_real_)) %>% setNames("grd_map")

##### env example #####
#set scenario
env_out_list <- list.files(glue::glue("{wd}environment/RCP45/2070"), pattern = ".tif")

pts          <- grd_out %>% st_as_sf() %>% st_centroid() %>% suppressWarnings() %>% st_geometry() %>%
                  st_transform(4326)

pred_pts     <- parallel::mclapply(env_out_list, function(x) {
  st <- read_stars(glue::glue("{wd}environment/RCP45/2070/{x}"))
  if(!grepl("bioc", x, fixed = T)) {
    st_extract(st, pts) %>% st_as_sf() %>% st_drop_geometry() %>% apply(., 1, median, na.rm = T)
  } else st_extract(st, pts) %>% st_as_sf() %>% st_drop_geometry()
}) %>% Reduce("cbind",.) %>% setNames(names(env_stars))

##### predict ####
pred_stars_glm <- grd_out
pred_stars_glm$grd_map[!is.na(pred_stars_glm$grd_map)]  <- predict(glm_mod, newdata = pred_pts)
pred_glm <- st_as_stars(pred_stars_glm$grd_map)
#write_stars(pred_glm, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_glm.tiff")

pred_stars_gam <- grd_out
pred_stars_gam$grd_map[!is.na(pred_stars_gam$grd_map)]  <- predict(gam_mod, newdata = pred_pts)
pred_gam <- st_as_stars(pred_stars_gam$grd_map)
#write_stars(pred_gam, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_gam.tiff")

pred_stars_mars <- grd_out
pred_stars_mars$grd_map[!is.na(pred_stars_mars$grd_map)]  <- predict(mars_mod, newdata = pred_pts)
pred_mars <- st_as_stars(pred_stars_mars$grd_map)
#write_stars(pred_mars, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_mars.tiff")

pred_stars_maxent <- grd_out
pred_stars_maxent$grd_map[!is.na(pred_stars_maxent$grd_map)]  <- predict(maxent_mod, pred_pts) 
pred_maxent <- st_as_stars(pred_stars_maxent$grd_map)
#write_stars(pred_maxent, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_maxent.tiff")





##### Ensemble model #####
#predictions
library(stars)
###### with raster #######
library(raster)
pred_glm2r = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_glm.tiff")
pred_gam2r = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_gam.tiff")
pred_mars2r = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_mars.tiff")
pred_maxent2r = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/pred_maxent.tiff")

plot(pred_glm2)
plot(pred_gam2r)
plot(pred_mars2r)
plot(pred_maxent2r)

models = stack(pred_glm2r, pred_gam2r, pred_mars2r, pred_maxent2r)
names(models) <- c("glm", "gam", "mars", "maxent")
mean_mod = mean(models)
plot(mean_mod)


###### with stars: 8 approaches which did not work #####

ggplot() +
  geom_stars(data = tt) +
  scale_fill_viridis_c(na.value = NA, name = "Probability") +
  geom_sf(data = map, fill = NA, linewidth = 0.6) +
  theme_void()

fls_path = "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/"

fls      <- list.files(fls_path, pattern = ".tiff$",  
                         recursive = TRUE, full.names = T)

tt <- read_stars(fls) %>% merge() %>% st_apply(., 1:2, median, future = T)
plot(tt)

write_stars(tt, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/prediction_test/ensemble_pred.tiff")



