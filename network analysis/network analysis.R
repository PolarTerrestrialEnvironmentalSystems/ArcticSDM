#network analysis

#### packages ####
library(stars)
library(tidyr)
library(terra)
library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

#### input ####
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM"
path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/Results"


#region
ecoreg <- st_read(glue::glue("{data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska

#get all result species names
flsl      <- list.dirs(path, 
                       recursive = FALSE, full.names = T)

speciesnames = as.data.frame(matrix(nrow = 1, ncol=1))

# Split the path by the "/" character
for (i in 1:length(flsl)){
path_parts <- strsplit(flsl[i], "/")[[1]]

# Extract the last part of the path, which contains the species name
speciesnames[i,1] <- path_parts[length(path_parts)]
}

colnames(speciesnames)[1] <- "species"


#for(species in speciesnames$species){
#for(sp in 1:nrow(speciesnames)){
  
  scen = 2
  year = 3
  
  #each species 1,3,6
  sp = 3
  species = speciesnames[sp,1]
  load(glue::glue("{data}/Results/{species}/Predictions/pred_stars.rda")) #future
  fut <- pred_stars[scen][,,,year]

  #crop future
  reg_t<- st_transform(reg, st_crs(fut))
  futc3 <- fut %>% st_crop(reg_t) %>% 
   setNames("p") %>%
    mutate(p = ifelse(p>0.5, 1, NA))
  # futca
  # st_as_sf(.)
  # futc1 <- subset(futc, futc$'2086' > 0.5)

 coordsxy3 <- as.data.frame(st_coordinates(futc3, dims = c("x", "y"))) %>%
   .[,1:2] %>%
 mutate(species = species)
 
 coordsxy6$species = species
 coordsxy6$p = 1
 
#unique(coordsxy)
nettab = rbind(coordsxy1, coordsxy6)

#for all together
nettab1 = nettab %>%
  
  #fut %>% st_crop(reg_t) %>%
  #as.data.frame(st_coordinates(., dims = c("x", "y"))) %>%
  
  #mutate(species = species) %>%
  mutate(x1 = round(X, -2)) %>%
  mutate(y1 = round(Y, -2)) %>%
  unite(., col = site, x1, y1, sep = "") %>% #combine to site
  .[,4:6] %>% #only species and site
  unique(.) #remove dublicates through rounding


# Transform the DataFrame to wide format
netwide <- nettab1 %>%
  pivot_wider(names_from = species, values_from = p,
              values_fn = list(p = sum)) 

summary(netwide)

netwide0 <- replace(netwide, is.na(netwide), 0)

#### create net ####
net_all <- graph_from_data_frame(netwide0, directed = FALSE)

#### clustering ####
cl <- cluster_louvain(net_all) 

#### output ####.
# check number of communities
length(cl) #95 communities detected

# check membership
membership(cl)

membercl <- data.frame(group = cl$membership, label = cl$names)