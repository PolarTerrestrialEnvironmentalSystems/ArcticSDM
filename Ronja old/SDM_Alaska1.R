library(sf)
sf_use_s2(FALSE)
library(stars)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(mgcv)
library(earth)

#data <- "/bioing/data/ArcticSDM_data/"

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/"

# ### Projection
# proj <- "+proj=laea +lon_0=-170 +lat_0=90"
# 
# ### Maps
# ecoreg   <- st_read(glue::glue("{data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
#   filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0) %>%
#   suppressWarnings() 
# 
# map_wrld <- st_read(glue::glue("{data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
#   st_intersection(ecoreg %>% st_union()) %>% dplyr::select(c('name', 'admin')) %>%
#   suppressWarnings()
# 
# 
# ### Subset
# map  <- map_wrld %>% filter(name=="Alaska")  %>% 
#   st_transform(proj)
# plot(map %>% dplyr::select(name))
# 
# ### Resolution
# res <- 5000


##############
## GBIF Occ ##
##############

# gbif_grid <- st_read(glue::glue("{map_data}occurance_grid.shp"))
# 
# ids <- which(c(gbif_grid %>% st_transform(st_crs(map)) %>%
#                  st_intersects(map, sparse = FALSE)))
# 
# ggplot() +
#   geom_sf(data = map) +
#   geom_sf(data = gbif_grid[ids,] %>% st_transform(st_crs(map)), fill = NA)
# 
# ### Read in data files
# gbif_list <- tibble(fls = list.files(glue::glue("{data}/GBIF/Plantea"))) %>%
#   mutate(ids = as.numeric(gsub("_gbif_plantea.csv", "", fls)))
# 
# ### Filter
# gbifTab    <- parallel::mclapply(which(gbif_list$ids%in%ids), function(x) {
#   read_csv(glue::glue("{data}/GBIF/Plantea/{gbif_list$fls[x]}"), progress = F, show_col_types = FALSE) %>%
#     dplyr::select(phylum, class, order, family, genus, species, scientificName, year, countryCode, 
#                   decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, basisOfRecord) %>%
#     filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species), 
#            is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 5000, year > 1970)
# }, mc.cores = parallel::detectCores()-5) %>% Reduce("rbind",.)
# 
# gbifTab_sf <- gbifTab %>% filter(!is.na(as.numeric(decimalLongitude))) %>%
#   st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
#   st_transform(st_crs(map))
# 
# ## Thinning
# res_thin <- 15000
# thin_map <- st_rasterize(map, st_as_stars(st_bbox(map), 
#                                           dx = res_thin, dy = res_thin, 
#                                           values = NA_real_)) %>% setNames("grd_map")
# 
# gbif_thin <- (gbifTab_sf %>% group_split(species))[sapply(gbifTab_sf %>% group_split(species), nrow)>50] %>%
#   parallel::mclapply(function(x) {
#     out <- x %>% dplyr::select(species) %>% st_transform(st_crs(thin_map)) %>%
#       mutate(cell = unlist(apply(st_intersects(., thin_map %>% st_as_sf(), sparse = FALSE), 1, function(y) ifelse(any(y), which(y), NA)))) %>%
#       filter(!is.na(cell)) 
#     if(nrow(out)>1) {
#       out %>% group_split(cell) %>% lapply(function(z) z %>% slice(sample(1:nrow(z), 1))) %>%
#         Reduce("rbind", .) %>% dplyr::select(cell, species) %>% suppressWarnings()
#     } else NULL
#   }, mc.cores = 8) %>% Reduce("rbind", .)


########################
## Environmental data ##
########################

# env_list  <- list.files(glue::glue("{wd}environment/calibration"), pattern = ".tiff")
# env_stars <- parallel::mclapply(env_list, function(x) {
#   st <- read_stars(glue::glue("{wd}environment/calibration/{x}"))
#   st_extract(st, gbif_thin %>% st_transform(st_crs(st))) %>% st_as_sf() %>% st_drop_geometry()
# }, mc.cores = length(env_list)) %>% Reduce("cbind", .) %>%
#   setNames(c(paste0("bio", 1:19), "prec", "tmax", "tmin"))
# 
# subsetList <- list(gbif = gbif_thin,
#                    env = env_stars)
#
#save(subsetList, file = glue::glue("{data}/GBIF/AlaskaList.rda"))
load(glue::glue("{data}/GBIF/AlaskaList.rda"))

##################
## Species Loop ##
##################

species <- subsetList$gbif %>% group_by(species) %>% st_drop_geometry() %>%
  summarise(sample = n()) %>% arrange(desc(sample))

# for(sp in species %>% pull(species)) {
sp <- (species %>% pull(species))[1]

occID <- subsetList$gbif %>% rownames_to_column(var = 'id') %>%
  filter(species==sp)

absID  <- subsetList$gbif %>% rownames_to_column(var = 'id') %>%
  filter(species!=sp, !duplicated(cell), !(cell%in%occID$cell))

modTab <- occID %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[occID$id,]) %>%
  bind_rows(absID %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[absID$id,]))

##### GLM #####
formula_glm <- paste0("p ~ ", paste(colnames(modTab)[-1], collapse = " + "))
glm_mod     <- glm(formula_glm, family = binomial(link = "logit"), data = modTab)
glm_pre     <- predict(glm_mod, type = "response")
plot(glm_pre)
glm_pre_s = st_as_stars(glm_pre)

##### GAM #####
#The formulae supplied to gam are exactly like those supplied to glm 
#except that smooth terms, s, te, ti and t2 can be added to the right hand side

library(mgcv)
formula_gam <- formula(paste0("p ~ ", paste0(paste0("s(", colnames(modTab)[-c(1)], collapse = ") + "), ")"))) 
formula_gam
gam_mod     <- gam(formula_gam, family = binomial(link = "logit"), data = modTab)  
gam_pre     <- predict(gam_mod, type = "response")
plot(gam_pre)
gam_pre_s = st_as_stars(gam_pre)

##### MARS #####
library(earth)
modTabc = na.omit(modTab)
mars_mod  <- earth(p ~ .,  data = modTabc, glm=list(family=binomial(link = "logit"))) 
    #needed to exclude NAs, needed to add the glm argument binomial
mars_pre <- predict(mars_mod, type = "response")
plot(mars_pre)
mars_pre_s = st_as_stars(mars_pre)

##### MaxEnt ###
## potentially add ENMevaluate
library(dismo)
args_list <- c('responsecurves=TRUE',
               'jackknife=TRUE',
               'pictures=TRUE',
               'autofeature=FALSE',
               'linear=TRUE',
               'quadratic=TRUE',
               'product=TRUE',
               'threshold=TRUE',
               'hinge=TRUE')

maxent_mod <- maxent(modTab[,-1], p = modTab$p, a=NULL, removeDuplicates=TRUE, nbg=10000,
                     args=args_list)
maxent_pre <- predict(maxent_mod, modTab, type = "response") #needed to input data
plot(maxent_pre)


# }
library(matrixStats)
library(dplyr)
all_pre <- as.data.frame(cbind(glm_pre, gam_pre, mars_pre, maxent_pre))
colnames(all_pre)[3] = "mars_pre"
all_pre_av2 <- all_pre  %>% mutate(mean = rowMeans(.[,1:4]))

###
 # mutate(median = rowMedians(.[,1:4]))
 # mutate(mean = mean(glm_pre, gam_pre, mars_pre, maxent_pre)) %>%
  #mutate(median = median(glm_pre, gam_pre, maxent_pre, maxent_pre))
 # median = as.matrix(.[grep('M\\d+', names(.))])) 


