#4step 

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

data <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/"
env  <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/"

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### Maps
load("//smb.isipd.dmawi.de/home/pd/slisovsk/Documents/ArcticSDM/GBIF_dataComp/gbif_regions.rda")

### Subset
plot(map %>% dplyr::select(region))

### Resolution
res <- 5000

########################
## GBIF Occ ############
########################

regMap    <- map %>% dplyr::select(region) %>% filter(region == "West")

gbif_thin <- tibble(files = list.files("//smb.isipd.dmawi.de/home/pd/slisovsk/Documents/ArcticSDM/GBIF_dataComp/", pattern = "thin")) %>%
  filter(grepl(regMap %>% pull(region) %>% unique() %>% as.character(), files)) %>%
  rownames_to_column(var = "row") %>% group_split(row) %>%
  lapply(function(x) {
    load(glue::glue("//smb.isipd.dmawi.de/home/pd/slisovsk/Documents/ArcticSDM/GBIF_dataComp/{x$files[1]}")); gbif_thin
  }) %>% Reduce("rbind", .)


########################
## Modern Environment ##
########################

wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
env_list  <- list.files(glue::glue("{wd}environment/CHELSA data 5 km/Modern/"), pattern = ".tif")

setwd(glue::glue("{wd}environment/CHELSA data 5 km/Modern/"))
getwd()
stars_env  <- read_stars(env_list) %>% setNames(sapply(strsplit(names(.), "_"), function(x) x[[1]]))

env_stars <- parallel::mclapply(env_list, function(x) {
  st <- read_stars(glue::glue("{wd}environment/CHELSA data 5 km/Modern/{x}"))
  st_extract(st, gbif_thin %>% st_transform(st_crs(st))) %>% st_as_sf() %>% st_drop_geometry()
}, mc.cores = 1) %>% Reduce("cbind", .) %>%
  setNames(c(paste0("bio", 1:19), "scd", "swe", "tasmax", "tasmin"))

subsetList <- list(gbif = gbif_thin,
                   env = env_stars)


########################
## Species Loop ########
########################

species <- subsetList$gbif %>% group_by(species) %>% st_drop_geometry() %>%
  summarise(sample = n()) %>% arrange(desc(sample))

# for(sp in species %>% pull(species)) {

group_by(species) %>% summarise(ncells = n()) %>%
  arrange(desc(ncells))


sp <- "Picea glauca"

occSpec <- gbif_thin %>% filter(species == sp)

occID <- subsetList$gbif %>% rownames_to_column(var = 'id') %>%
  filter(species==sp)

absID  <- subsetList$gbif %>% rownames_to_column(var = 'id') %>%
  filter(species!=sp, !duplicated(cell), !(cell%in%occID$cell))

modTab <- occID %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[occID$id,]) %>%
  bind_rows(absID %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[absID$id,]))

##### Step 1: Specifying geographical extent  #####

##################
## Species Loop ##
##################

presence_sf <- st_as_sf(occID, coords = c("geometry"), crs = 4326) 
absence_sf <- st_as_sf(absID, coords = c("geometry"), crs = 4326)

####################################
###### create buffers 200km ######
presence_buffer200 <- st_buffer(presence_sf, dist = 200000) 

intersection200abs <- st_intersection(absence_sf, presence_buffer200)
# select the absence data from the buffer

buffer200abs <- intersection200abs[!duplicated(st_coordinates(intersection200abs)), ]
#only selected the pressumed abs inside the buffer around pres

#create data tab
modTab_buf200 <- occID %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[occID$id,]) %>%
  bind_rows(buffer200abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer200abs$id,]))


rpts = na.omit(modTab_buf200)

rpts_pres = subset(rpts, rpts$p ==1)


#### step2 ####

library(e1071)
mod <- svm(rpts[,-1], rpts[,1], type = "one-classification")
backg <- rpts %>% filter(p==0 & predict(mod)==0) %>% dplyr::select(-p)

modelTab <- backg %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) %>% bind_rows(rpts_pres) %>% na.omit()


#### step3 ####
km <- kmeans(backg, centers = nrow(rpts_pres))
modelTab <- km$centers %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) %>% bind_rows(rpts_pres) %>% na.omit()


#### glm ####
glm_formula <- paste0("p ~ ", paste(names(modelTab[,-1]), collapse = " + "))
glmMod <- glm(glm_formula, data = modelTab, family = binomial(link = "logit"))

pred_glm <- predict(glmMod,  type = "response")

plot(pred_glm, cex = 0.5, pch = 16, col = ifelse(modelTab$p==1, "red", "blue"))
plot(data.frame(g = as.factor(modelTab$p), p = pred_glm))

tibble(g = modelTab$p, p = pred_glm) %>% group_by(g) %>% summarise(median = median(p))


pROC::auc(pROC::roc(modelTab %>% pull(p), pred_glm))

#}

