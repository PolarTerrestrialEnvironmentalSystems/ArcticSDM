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
sp <- "Picea glauca"



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


##### GAM #####
#The formulae supplied to gam are exactly like those supplied to glm 
#except that smooth terms, s, te, ti and t2 can be added to the right hand side

library(mgcv)
formula_gam <- formula(paste0("p ~ ", paste0(paste0("s(", colnames(modTab)[-c(1)], collapse = ") + "), ")"))) 
formula_gam
gam_mod     <- gam(formula_gam, family = binomial(link = "logit"), data = modTab)  
gam_pre     <- predict(gam_mod, type = "response")
plot(gam_pre)


##### MARS #####
library(earth)
modTabc = na.omit(modTab)
mars_mod  <- earth(p ~ .,  data = modTabc, glm=list(family=binomial(link = "logit"))) 
    #needed to exclude NAs, needed to add the glm argument binomial
mars_pre <- predict(mars_mod, type = "response")
plot(mars_pre)


##### MaxEnt ####
## potentially add ENMevaluate
library(dismo)
args_list <- c('responsecurves=TRUE',
               'jackknife=TRUE',
               'pictures=TRUE',
               'autofeature=TRUE')
               
maxent_mod <- maxent(modTab[,-1], p = modTab$p, a=NULL, removeDuplicates=TRUE, nbg=0,
                     args=args_list)
maxent_pre <- predict(maxent_mod, modTab, type = "response") #needed to input data
plot(maxent_pre)


# }
library(matrixStats)
library(dplyr)
all_pre <- as.data.frame(cbind(glm_pre, gam_pre, mars_pre, maxent_pre))
colnames(all_pre)[3] = "mars_pre"
all_pre_av <- all_pre  %>% mutate(mean = rowMeans(.[,1:4]))

#### auc ####
library(pROC)
auc1 <- auc(roc(modTabc$p, glm_pre))
auc2 <- auc(roc(modTabc$p, gam_pre))
auc3 <- auc(roc(modTabc$p, mars_pre))
auc4 <- auc(roc(modTabc$p, maxent_pre))

our_mod = c(auc1, auc2, auc3, auc4)
aucs = as.data.frame(our_mod)
aucs$mod = names(all_pre)            
aucs$ssdm = c(0.8128, 0.8698, 0.8353, 0.8624)
                                                               
boxplot(glm_mod$y, glm_mod$fitted.values)
t.test(glm_mod$y, glm_mod$fitted.values)
plot(glm_mod$fitted.values, glm_mod$y)
plot(glm_pre, glm_mod$y)

boxplot(modTabc$p, glm_pre)
tglm = t.test(modTabc$p, glm_pre)

boxplot(modTabc$p, gam_pre)
tgam = t.test(modTabc$p, gam_pre)

boxplot(modTabc$p, mars_pre)
tmars = t.test(modTabc$p, mars_pre)

boxplot(modTabc$p, maxent_pre)
tmaxent = t.test(modTabc$p, maxent_pre)

aucs$our_p = c(tglm$p.value, tgam$p.value, tmars$parameter, tmaxent$p.value)

AUC <- as.data.frame(aucs[with(aucs,
                               order(aucs[,"our_mod"], decreasing=TRUE)), ])
tmaxent

#### 3steps method ####
##### buffering #####
presence_sf <- st_as_sf(occID, coords = c("geometry"), crs = 4326) 
absence_sf <- st_as_sf(absID, coords = c("geometry"), crs = 4326)

####################################

#### short ####

#there is a bracket error, and I can see why, but I want the bracket for the
#lapply loop to go around everything

buf_list <- list(1, 5, 10, 20, 50, 100, 200, 500, 5000)

lapply(1:length(buf_list), function(f) 
  
bufferabs <- st_buffer(presence_sf, dist = buf_list[f]) %>%
  st_intersection(absence_sf, .) %>%
  .[!duplicated(st_coordinates(.)),] #intersection5000pres[!duplicated(st_coordinates(intersection5000pres)), ]

modTab_buf <- occID %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[occID$id,]) %>%
  bind_rows(bufferabs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[bufferabs$id,])))

library(FactoMineR)

results <- list()
results[[f]] <- 
  PCA(na.omit(modTab_buf), ncp=5) %>% #pca
  
  as.data.frame(res.pca5000$var$contrib) %>% #subset 
  .[with(., order(.[,"Dim.1"], decreasing=TRUE)),] %>% #sort
  
  as.data.frame(.[1:3, 1:2]) %>% #subset
  mutate(var = rownames(.), #add infos
         buf = buf_list[f]))

###### 1 #####

###### create buffers ######
presence_buffer1 <- st_buffer(presence_sf, dist = 1000) 

intersection1abs <- st_intersection(absence_sf, presence_buffer1)
# select the absence data from the buffer

buffer1abs <- intersection1abs[!duplicated(st_coordinates(intersection1abs)), ]


#create data tab
modTab_buf1a <- occID %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[occID$id,]) %>%
  bind_rows(buffer1abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer1abs$id,]))

###### pca  ######
modTab_buf1 = na.omit(modTab_buf1a)


library(FactoMineR)
res.pca1 <- PCA(modTab_buf1, ncp=5) # Compute PCA


###### variable selection #####
var_cont01 <- as.data.frame(res.pca1$var$contrib)
var_cont1 <- as.data.frame(var_cont01[with(var_cont01,
                                           order(var_cont01[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub1 = as.data.frame(var_cont1[1:3, 1:2])
var_cont_sub1$var = rownames(var_cont_sub1)
var_cont_sub1$buf = 1

#################

###### 5 #####

###### create buffers ######
presence_buffer5 <- st_buffer(presence_sf, dist = 5000) 

intersection5abs <- st_intersection(absence_sf, presence_buffer5)
# select the absence data from the buffer

buffer5abs <- intersection5abs[!duplicated(st_coordinates(intersection5abs)), ]


intersection5pres <- st_intersection(presence_sf, presence_buffer5)
# select the presence data from the buffer


buffer5pres <- intersection5pres[!duplicated(st_coordinates(intersection5pres)), ]

#create data tab
modTab_buf5 <- buffer5pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer5pres$id,]) %>%
  bind_rows(buffer5abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer5abs$id,]))

###### pca  ######
modTab_buf5 = na.omit(modTab_buf5)


library(FactoMineR)
res.pca5 <- PCA(modTab_buf5, ncp=5) # Compute PCA


###### variable selection #####
var_cont05 <- as.data.frame(res.pca5$var$contrib)
var_cont5 <- as.data.frame(var_cont05[with(var_cont05,
                                           order(var_cont05[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub5 = as.data.frame(var_cont5[1:3, 1:2])
var_cont_sub5$var = rownames(var_cont_sub5)
var_cont_sub5$buf = 5

###### 10 #####

###### create buffers ######
presence_buffer10 <- st_buffer(presence_sf, dist = 10000) 

intersection10abs <- st_intersection(absence_sf, presence_buffer10)
# select the absence data from the buffer

buffer10abs <- intersection10abs[!duplicated(st_coordinates(intersection10abs)), ]


intersection10pres <- st_intersection(presence_sf, presence_buffer10)
# select the presence data from the buffer


buffer10pres <- intersection10pres[!duplicated(st_coordinates(intersection10pres)), ]

#create data tab
modTab_buf10 <- buffer10pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer10pres$id,]) %>%
  bind_rows(buffer10abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer10abs$id,]))

###### pca  ######
modTab_buf10 = na.omit(modTab_buf10)
10

library(FactoMineR)
res.pca10 <- PCA(modTab_buf10, ncp=5) # Compute PCA


###### variable selection #####
var_cont010 <- as.data.frame(res.pca10$var$contrib)
var_cont10 <- as.data.frame(var_cont010[with(var_cont010,
                                           order(var_cont010[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub10 = as.data.frame(var_cont10[1:3, 1:2])
var_cont_sub10$var = rownames(var_cont_sub10)
var_cont_sub10$buf = 10

###### 20 ######

###### create buffers ######
presence_buffer20 <- st_buffer(presence_sf, dist = 20000) 
#=20 km? or 20000*5km because 5000 is resolution

intersection20abs <- st_intersection(absence_sf, presence_buffer20)
# select the absence data from the buffer

buffer20abs <- intersection20abs[!duplicated(st_coordinates(intersection20abs)), ]


intersection20pres <- st_intersection(presence_sf, presence_buffer20)
# select the presence data from the buffer


buffer20pres <- intersection20pres[!duplicated(st_coordinates(intersection20pres)), ]

#create data tab
modTab_buf20 <- buffer20pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer20pres$id,]) %>%
  bind_rows(buffer20abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer20abs$id,]))

###### pca  ######
modTab_buf20 = na.omit(modTab_buf20)
#res.pca20 <- prcomp(modTab_buf, scale = TRUE) # Compute PCA
#res.pca2 <- princomp(modTab_buf) # Compute PCA

library(FactoMineR)
res.pca20 <- PCA(modTab_buf20, ncp=5) # Compute PCA
#res.pca20$var$contrib
#note: there is 49% and 10% of variables in dim5

###### variable selection #####
var_cont020 <- as.data.frame(res.pca20$var$contrib)
var_cont20 <- as.data.frame(var_cont020[with(var_cont020,
                        order(var_cont020[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub20 = as.data.frame(var_cont20[1:3, 1:2])
var_cont_sub20$var = rownames(var_cont_sub20)
var_cont_sub20$buf = 20


############################
###### 50 #####

###### create buffers ######
presence_buffer50 <- st_buffer(presence_sf, dist = 50000) 

intersection50abs <- st_intersection(absence_sf, presence_buffer50)
# select the absence data from the buffer

buffer50abs <- intersection50abs[!duplicated(st_coordinates(intersection50abs)), ]


intersection50pres <- st_intersection(presence_sf, presence_buffer50)
# select the presence data from the buffer


buffer50pres <- intersection50pres[!duplicated(st_coordinates(intersection50pres)), ]

#create data tab
modTab_buf50 <- buffer50pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer50pres$id,]) %>%
  bind_rows(buffer50abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer50abs$id,]))

###### pca  ######
modTab_buf50 = na.omit(modTab_buf50)


library(FactoMineR)
res.pca50 <- PCA(modTab_buf50, ncp=5) # Compute PCA


###### variable selection #####
var_cont050 <- as.data.frame(res.pca50$var$contrib)
var_cont50 <- as.data.frame(var_cont050[with(var_cont050,
                                           order(var_cont050[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub50 = as.data.frame(var_cont50[1:3, 1:2])
var_cont_sub50$var = rownames(var_cont_sub50)
var_cont_sub50$buf = 50

################################################

###### 100 #####

###### create buffers ######
presence_buffer100 <- st_buffer(presence_sf, dist = 100000) 

intersection100abs <- st_intersection(absence_sf, presence_buffer100)
# select the absence data from the buffer

buffer100abs <- intersection100abs[!duplicated(st_coordinates(intersection100abs)), ]


intersection100pres <- st_intersection(presence_sf, presence_buffer100)
# select the presence data from the buffer


buffer100pres <- intersection100pres[!duplicated(st_coordinates(intersection100pres)), ]

#create data tab
modTab_buf100 <- buffer100pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer100pres$id,]) %>%
  bind_rows(buffer100abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer100abs$id,]))

###### pca  ######
modTab_buf100 = na.omit(modTab_buf100)
10

library(FactoMineR)
res.pca100 <- PCA(modTab_buf100, ncp=5) # Compute PCA


###### variable selection #####
var_cont0100 <- as.data.frame(res.pca100$var$contrib)
var_cont100 <- as.data.frame(var_cont0100[with(var_cont0100,
                                           order(var_cont0100[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub100 = as.data.frame(var_cont100[1:3, 1:2])
var_cont_sub100$var = rownames(var_cont_sub100)
var_cont_sub100$buf = 100

#########################

###### 200 #####

###### create buffers ######
presence_buffer200 <- st_buffer(presence_sf, dist = 200000) 

intersection200abs <- st_intersection(absence_sf, presence_buffer200)
# select the absence data from the buffer

buffer200abs <- intersection200abs[!duplicated(st_coordinates(intersection200abs)), ]


intersection200pres <- st_intersection(presence_sf, presence_buffer200)

###############################
# select the presence data from the buffer


buffer200pres <- intersection200pres[!duplicated(st_coordinates(intersection200pres)), ]

#create data tab
modTab_buf200 <- buffer200pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer200pres$id,]) %>%
  bind_rows(buffer200abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer200abs$id,]))

###### pca  ######
modTab_buf200 = na.omit(modTab_buf200)
10

library(FactoMineR)
res.pca200 <- PCA(modTab_buf200, ncp=5) # Compute PCA


###### variable selection #####
var_cont0200 <- as.data.frame(res.pca200$var$contrib)
var_cont200 <- as.data.frame(var_cont0200[with(var_cont0200,
                                           order(var_cont0200[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub200 = as.data.frame(var_cont200[1:3, 1:2])
var_cont_sub200$var = rownames(var_cont_sub200)
var_cont_sub200$buf = 200

##############################

###### 300 #####

###### create buffers ######
presence_buffer300 <- st_buffer(presence_sf, dist = 300000) 

intersection300abs <- st_intersection(absence_sf, presence_buffer300)
# select the absence data from the buffer

buffer300abs <- intersection300abs[!duplicated(st_coordinates(intersection300abs)), ]


intersection300pres <- st_intersection(presence_sf, presence_buffer300)

###############################
# select the presence data from the buffer


buffer300pres <- intersection300pres[!duplicated(st_coordinates(intersection300pres)), ]

#create data tab
modTab_buf300 <- buffer300pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer300pres$id,]) %>%
  bind_rows(buffer300abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer300abs$id,]))

###### pca  ######
modTab_buf300 = na.omit(modTab_buf300)
10

library(FactoMineR)
res.pca300 <- PCA(modTab_buf300, ncp=5) # Compute PCA


###### variable selection #####
var_cont0300 <- as.data.frame(res.pca300$var$contrib)
var_cont300 <- as.data.frame(var_cont0300[with(var_cont0300,
                                               order(var_cont0300[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub300 = as.data.frame(var_cont300[1:3, 1:2])
var_cont_sub300$var = rownames(var_cont_sub300)
var_cont_sub300$buf = 300

##############################

###### 400 #####

###### create buffers ######
presence_buffer400 <- st_buffer(presence_sf, dist = 400000) 

intersection400abs <- st_intersection(absence_sf, presence_buffer400)
# select the absence data from the buffer

buffer400abs <- intersection400abs[!duplicated(st_coordinates(intersection400abs)), ]


intersection400pres <- st_intersection(presence_sf, presence_buffer400)

###############################
# select the presence data from the buffer


buffer400pres <- intersection400pres[!duplicated(st_coordinates(intersection400pres)), ]

#create data tab
modTab_buf400 <- buffer400pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer400pres$id,]) %>%
  bind_rows(buffer400abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer400abs$id,]))

###### pca  ######
modTab_buf400 = na.omit(modTab_buf400)
10

library(FactoMineR)
res.pca400 <- PCA(modTab_buf400, ncp=5) # Compute PCA


###### variable selection #####
var_cont0400 <- as.data.frame(res.pca400$var$contrib)
var_cont400 <- as.data.frame(var_cont0400[with(var_cont0400,
                                               order(var_cont0400[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub400 = as.data.frame(var_cont400[1:3, 1:2])
var_cont_sub400$var = rownames(var_cont_sub400)
var_cont_sub400$buf = 400

##############################
###### 500 #####

###### create buffers ######
presence_buffer500 <- st_buffer(presence_sf, dist = 500000) 

intersection500abs <- st_intersection(absence_sf, presence_buffer500)
# select the absence data from the buffer

buffer500abs <- intersection500abs[!duplicated(st_coordinates(intersection500abs)), ]


intersection500pres <- st_intersection(presence_sf, presence_buffer500)
# select the presence data from the buffer


buffer500pres <- intersection500pres[!duplicated(st_coordinates(intersection500pres)), ]

#create data tab
modTab_buf500 <- buffer500pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer500pres$id,]) %>%
  bind_rows(buffer500abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer500abs$id,]))

###### pca  ######
modTab_buf500 = na.omit(modTab_buf500)
10

library(FactoMineR)
res.pca500 <- PCA(modTab_buf500, ncp=5) # Compute PCA


###### variable selection #####
var_cont0500 <- as.data.frame(res.pca500$var$contrib)
var_cont500 <- as.data.frame(var_cont0500[with(var_cont0500,
                                           order(var_cont0500[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub500 = as.data.frame(var_cont500[1:3, 1:2])
var_cont_sub500$var = rownames(var_cont_sub500)
var_cont_sub500$buf = 500


###### 5000 #####

###### create buffers ######
presence_buffer5000 <- st_buffer(presence_sf, dist = 5000000) 

intersection5000abs <- st_intersection(absence_sf, presence_buffer5000)
# select the absence data from the buffer

buffer5000abs <- intersection5000abs[!duplicated(st_coordinates(intersection5000abs)), ]


intersection5000pres <- st_intersection(presence_sf, presence_buffer5000)
# select the presence data from the buffer


buffer5000pres <- intersection5000pres[!duplicated(st_coordinates(intersection5000pres)), ]

#create data tab
modTab_buf5000 <- buffer5000pres %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[buffer5000pres$id,]) %>%
  bind_rows(buffer5000abs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[buffer5000abs$id,]))

###### pca  ######
modTab_buf5000 = na.omit(modTab_buf5000)


library(FactoMineR)
res.pca5000 <- PCA(modTab_buf5000, ncp=5) # Compute PCA


###### variable selection #####
var_cont05000 <- as.data.frame(res.pca5000$var$contrib)
var_cont5000 <- as.data.frame(var_cont05000[with(var_cont05000,
                                            order(var_cont05000[,"Dim.1"], decreasing=TRUE)), ])
var_cont_sub5000 = as.data.frame(var_cont5000[1:3, 1:2])
var_cont_sub5000$var = rownames(var_cont_sub5000)
var_cont_sub5000$buf = 5000

##### results #####
results <- rbind(var_cont_sub1, var_cont_sub5, var_cont_sub20, var_cont_sub50,
                 var_cont_sub100, var_cont_sub200,  var_cont_sub300, var_cont_sub400,
                 var_cont_sub500, var_cont_sub5000)

results.prec = subset(results, results$var =="prec")
plot(results.prec$buf, results.prec$Dim.1, type = "o")


results.bio14 = subset(results, results$var =="bio14")
plot(results.prec$buf, results.prec$Dim.1, type = "o")
lines(results.bio14$buf, results.bio14$Dim.1, type = "o", col="red")

best_buf <- results.prec$buf[which.max(results.prec$Dim.1)]


#### short ####

#verschachtelt base
#st_intersection(absence_sf, (st_buffer(presence_sf, dist = 1000)))
#buffer1abs <- st_intersection(absence_sf, (st_buffer(presence_sf, dist = 1000)))[!duplicated(st_coordinates(st_intersection(absence_sf, (st_buffer(presence_sf, dist = 1000))))), ]

#there is a bracket error, and I can see why, but I want the bracket for the
#lapply loop to go around everything

buf_list <- list(1, 5, 10, 20, 50, 100, 200, 500, 5000)

lapply(1:length(buf_list), function(f) 
  
bufferabs <- st_buffer(presence_sf, dist = buf_list[f]) %>%
  st_intersection(absence_sf, .) %>%
  .[!duplicated(st_coordinates(.)),] #intersection5000pres[!duplicated(st_coordinates(intersection5000pres)), ]

modTab_buf <- occID %>% mutate(p = 1) %>% dplyr::select(p) %>% st_drop_geometry() %>%
  bind_cols(subsetList$env[occID$id,]) %>%
  bind_rows(bufferabs %>% mutate(p = 0) %>% dplyr::select(p) %>% st_drop_geometry() %>%
              bind_cols(subsetList$env[bufferabs$id,])))

library(FactoMineR)

results <- list()
results[[f]] <- 
  PCA(na.omit(modTab_buf), ncp=5) %>% #pca
  
  as.data.frame(res.pca5000$var$contrib) %>% #subset 
  .[with(., order(.[,"Dim.1"], decreasing=TRUE)),] %>% #sort
  
  as.data.frame(.[1:3, 1:2]) %>% #subset
  mutate(var = rownames(.), #add infos
         buf = buf_list[f])) 





