library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)


# wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
  wd <- "/Volumes/projects/bioing/data/ArcticSDM/"
# wd <- '/bioing/data/ArcticSDM/'
  map_data <- "~/Documents/ArtcticSDM_data/"
# map_data <- '/bioing/data/ArcticSDM/data/'

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"
  
### Maps
ecoreg   <- st_read(glue::glue("{map_data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

map_wrld <- st_read(glue::glue("{map_data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
  st_intersection(ecoreg %>% st_union()) %>% dplyr::select(c('name', 'admin'))


### Subset
map <- map_wrld %>% filter(name=="Alaska")  %>% 
  st_transform(proj)

### Resolution
res <- 5000
grd_map <- st_rasterize(map, st_as_stars(st_bbox(map), 
             dx = res, dy = res, 
             values = NA_real_)) %>% setNames("grd_map")


##############
## GBIF Occ ##
##############

gbif_grid <- st_read(glue::glue("{map_data}occurance_grid.shp"))

ids <- which(c(gbif_grid %>% st_transform(st_crs(map)) %>%
  st_intersects(map, sparse = FALSE)))

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = gbif_grid[ids,] %>% st_transform(st_crs(map)), fill = NA)
 
### Read in data files
gbif_list <- tibble(fls = list.files(glue::glue("{wd}/data/GBIF/Plantea"))) %>%
  mutate(ids = as.numeric(gsub("_gbif_plantea.csv", "", fls)))

### Filter
gbifTab    <- parallel::mclapply(which(gbif_list$ids%in%ids), function(x) {
  read_csv(glue::glue("{wd}Data/GBIF/Plantea/{gbif_list$fls[x]}"), progress = F, show_col_types = FALSE) %>%
    dplyr::select(kingdom, phylum, class, order, family, genus, species, scientificName, countryCode, 
                  decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters) %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species)) 
  }, mc.cores = 8) %>% Reduce("rbind",.) %>%
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
  st_transform(st_crs(map))

## Thinning
res_thin <- 15000
thin_map <- st_rasterize(map, st_as_stars(st_bbox(map), 
                                         dx = res_thin, dy = res_thin, 
                                         values = NA_real_)) %>% setNames("grd_map")

gbif_thin <- (gbifTab %>% group_split(species))[sapply(gbifTab %>% group_split(species), nrow)>50] %>%
  parallel::mclapply(function(x) {
    st_rasterize(x %>% dplyr::select(species), thin_map) %>% st_as_sf(na.rm = FALSE) %>%
      mutate(poly = st_as_sf(thin_map, na.rm = FALSE) %>% pull(grd_map)) %>%
      st_centroid() %>%
      rownames_to_column(var = "id") %>%
      filter(!is.na(poly) & !is.na(ID) & ID > 1) %>%
      mutate(species = x$species[1]) %>%
      dplyr::select(id, species) %>% suppressWarnings()
    }, mc.cores = 8) %>% Reduce("rbind", .)


########################
## Environmental data ##
########################

tt <- gbif_thin %>% dplyr::select(id) %>% filter(!duplicated(id))

env_list  <- list.files(glue::glue("{wd}environment/calibration"), pattern = ".tiff")
env_stars <- lapply(env_list, function(x) {
  st <- read_stars(glue::glue("{wd}environment/calibration/{x}"))
  st_extract(st, tt %>% st_transform(st_crs(st))) %>% st_as_sf() %>% st_drop_geometry()
  }) %>% Reduce("cbind", .)


##################
## Species Loop ##
##################

sp <- gbif_thin %>% st_drop_geometry() %>% group_by(species) %>% summarise(sample = n()) %>%
  arrange(desc(sample))

s <- sp %>% slice(1) %>% pull(species)

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
  

### SMD ###

## GLM
formula <- paste0("p ~ ", paste(colnames(mod_tab)[-1], collapse = " + "))
glm_mod <- glm(formula, family = binomial(link = "logit"), data = mod_tab)


##### other models
#####
##### TO DO!
#####


### Predictions
res_pred <- 2500
grd_out  <- st_rasterize(map, st_as_stars(st_bbox(map), 
                                         dx = res_pred, dy = res_pred, 
                                         values = NA_real_)) %>% setNames("grd_map")

## env example
env_out_list <- list.files(glue::glue("{wd}environment/RCP45/2070"), pattern = ".tif")

pts          <- grd_out %>% st_as_sf() %>% st_centroid() %>% suppressWarnings() %>% st_geometry() %>%
                  st_transform(4326)

pred_pts     <- parallel::mclapply(env_out_list, function(x) {
  st <- read_stars(glue::glue("{wd}environment/RCP45/2070/{x}"))
  if(!grepl("bioc", x, fixed = T)) {
    st_extract(st, pts) %>% st_as_sf() %>% st_drop_geometry() %>% apply(., 1, median, na.rm = T)
  } else st_extract(st, pts) %>% st_as_sf() %>% st_drop_geometry()
}) %>% Reduce("cbind",.) %>% setNames(names(env_stars))

pred_stars <- grd_out
pred_stars$grd_map[!is.na(pred_stars$grd_map)]  <- predict(glm_mod, newdata = pred_pts)

ggplot() +
  geom_stars(data = pred_stars) +
  scale_fill_viridis_c(na.value = NA, name = "Probability") +
  geom_sf(data = map, fill = NA, linewidth = 0.6) +
  theme_void()
  
