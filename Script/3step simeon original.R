#simeon 3step original

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
#load("//smb.isipd.dmawi.de/projects/bioing/data/GBIF_dataComp/gbif_regions.rda")
load("//smb.isipd.dmawi.de/home/pd/slisovsk/Documents/ArcticSDM/GBIF_dataComp/gbif_regions.rda")
### Subset
plot(map %>% dplyr::select(region))

### Resolution
res <- 5000

########################
## Modern Environment ##
########################

env_list   <- list.files(glue::glue("{env}/CHELSA data 5 km/Modern/"), pattern = ".tif", full.names = T)
stars_env  <- read_stars(env_list) %>% setNames(sapply(strsplit(names(.), "_"), function(x) x[[1]]))

########################
## GBIF Occ ############
########################

regMap    <- map %>% dplyr::select(region)  %>% filter(region == "West")

gbif_thin <- tibble(files = list.files("//smb.isipd.dmawi.de/home/pd/slisovsk/Documents/ArcticSDM/GBIF_dataComp/", pattern = "thin")) %>%
  filter(grepl(regMap %>% pull(region) %>% unique() %>% as.character(), files)) %>%
  rownames_to_column(var = "row") %>% group_split(row) %>%
  lapply(function(x) {
    load(glue::glue("//smb.isipd.dmawi.de/home/pd/slisovsk/Documents/ArcticSDM/GBIF_dataComp/{x$files[1]}")); gbif_thin
  }) %>% Reduce("rbind", .)


########################
## Species Loop ########
########################

species <- gbif_thin %>% st_drop_geometry() %>%
  group_by(species) %>% summarise(ncells = n()) %>%
  arrange(desc(ncells))

# for(sp in species %>% pull(species)) {

# sp <- "Picea glauca"

occSpec <- gbif_thin %>% filter(species == sp)

### EnvVariables
occTab <- st_extract(stars_env, occSpec) %>% st_drop_geometry() %>%
  mutate(p = 1, .before = names(.)[1])

### Pseudo absence
## Step 1
radius <- 50000
sample <- 5000

# pcaRadii <- parallel::mclapply(radius, function(x) {
#  
#     bfr <- occSpec %>% st_buffer(x) %>% st_union()
#     rpts <- st_sample(bfr, sample, type = "random") %>%
#       st_extract(stars_env, .) %>% st_drop_geometry() %>% na.omit()
#    
#     pca <- vegan::rda(rpts)
#     (100*vegan::scores(pca, display = "sp", scaling = 0)[,1]^2) %>% as_tibble() %>%
#       mutate(var = rownames(vegan::scores(pca)$species), radius = x)
#
# }, mc.cores = length(radius)) %>% Reduce("rbind", .)

bfr <- occSpec %>% st_buffer(radius) %>% st_union()

ggplot() +
  geom_sf(data = occSpec %>% st_geometry()) +
  geom_sf(data = bfr, fill = "transparent") +
  geom_sf(data = regMap %>% st_geometry(), fill = "transparent")

rpts <- st_sample(bfr, sample, type = "random") %>%
  st_extract(stars_env, .) %>% st_drop_geometry() %>% mutate(p = 0, .before = names(.)[1]) %>%
  bind_rows(occTab) %>% na.omit()

library(e1071)
mod <- svm(rpts[,-1], rpts[,1], type = "one-classification")
backg <- rpts %>% filter(p==0 & predict(mod)==0) %>% dplyr::select(-p)

km <- kmeans(backg, centers = nrow(occTab))
modelTab <- km$centers %>% as_tibble() %>% mutate(p = 0, .before = names(.)[1]) %>%
  bind_rows(occTab) %>% na.omit()


### glm
glm_formula <- paste0("p ~ ", paste(names(modelTab[,-1]), collapse = " + "))
glmMod <- glm(glm_formula, data = modelTab, family = binomial(link = "logit"))

pred_glm <- predict(glmMod, type = "response")

plot(pred_glm, cex = 0.5, pch = 16, col = ifelse(modelTab$p==1, "red", "blue"))
plot(data.frame(g = as.factor(modelTab$p), p = pred_glm))

tibble(g = modelTab$p, p = pred_glm) %>% group_by(g) %>% summarise(median = median(p))

newdata  <- stars_env[regMap %>% st_transform(st_crs(stars_env)),,,] %>% merge() %>% st_as_sf()
pred_mod <- newdata %>% mutate(pred = predict(glmMod, newdata = newdata %>% st_drop_geometry(), type = "response")) %>% dplyr::select(pred) %>%
  st_rasterize(., stars_env[regMap %>% st_transform(st_crs(stars_env)),,,][1])

ggplot() +
  geom_stars(data = pred_mod) +
  scale_fill_viridis_c() +
  geom_sf(data = occSpec$geometry, col = "magenta")

pROC::auc(pROC::roc(modelTab %>% pull(p), pred_glm))

# }
