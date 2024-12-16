library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)
library(flexsdm)
library(e1071)
library(pROC)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(dplyr)
library(terra)

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### Maps
ecoreg   <- st_read("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/Ecoregions/tnc_terr_ecoregions.shp") %>%
filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)


# Define the polar projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"


# Get world map data and transform it to WGS84
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, crs = 4326)


#### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[2] #tundra
taiga_ne <- hcl.colors(6, "Purples")[2] #america
taiga_eu <- hcl.colors(6, "Blues 3")[2] #europe

realms <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "white") +  # Background map
  geom_sf(data = ecoreg, aes(color = WWF_REALM2), show.legend = TRUE) +  # Map ecoregions by WWF_REALM2
  scale_color_manual(values = c(taiga_ne, taiga_eu))+  # Replace with your realm names and desired colors
                                
  coord_sf(crs = proj, 
           xlim = c(-5000000, 5000000),  # Adjust these values based on your projection
           ylim = c(-5000000, 4000000),  # Adjust these values based on your projection
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Ecoregions by WWF Realms",
       color = "realms")  # Add title and legend label


biomes <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "white") +  # Background map
  geom_sf(data = ecoreg, aes(color = WWF_MHTNAM), show.legend = TRUE) +  # Map ecoregions by WWF_REALM2
  scale_color_manual(values = c("blue", tundra))+  # Use a color scale that works well for categorical data
  coord_sf(crs = proj, 
           xlim = c(-5000000, 5000000),  # Adjust these values based on your projection
           ylim = c(-5000000, 4000000),  # Adjust these values based on your projection
           expand = FALSE) +
  theme_minimal() +
  labs(title = "Ecoregions by WWF Biomes",
       color = "biomes")  # Add title and legend label


library(cowplot)
combined_plot_res <- plot_grid(realms, biomes , ncol = 2)
print(combined_plot_res)

print(biomes)
print(realms)

nearctic_taiga <- ecoreg %>% filter(WWF_REALM2 == "Nearctic") %>% filter(WWF_MHTNAM == "Boreal Forests/Taiga")
palearctic_taiga <- ecoreg %>% filter(WWF_REALM2 == "Palearctic") %>% filter(WWF_MHTNAM == "Boreal Forests/Taiga")  
tundra <- ecoreg %>% filter(WWF_MHTNAM == "Tundra")  

setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering")
write.csv(nearctic_taiga,"nearctic_taiga.csv",  row.names = F)
write.csv(palearctic_taiga,"palearctic_taiga.csv",  row.names = F)
write.csv(tundra,"tundra.csv",  row.names = F)

setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/timeline_viridis")
membercl <- read.csv("membercl1.csv")


##### join ####

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"
load(glue::glue("{data}/sigma1Array.rda")) #array
load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid

#### taiga_ne ####

# Ensure nearctic_taiga is an sf object

nearctic_taiga <- ecoreg %>% filter(WWF_REALM2 == "Nearctic") %>% filter(WWF_MHTNAM == "Boreal Forests/Taiga")

com1 <- spArrayMeta[[2]] %>%
  mutate(sp = apply(spArray[,,1,1], 2, sum)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))


#same crs
nearctic_taiga_sf <- st_as_sf(nearctic_taiga)
rast_com1 <- st_as_sf(com1)
crs_nearctic <- st_crs(nearctic_taiga_sf)
crs_rast_com1 <- st_crs(rast_com1)

if (!st_crs(nearctic_taiga_sf) == st_crs(rast_com1)) {
  # Transform nearctic_taiga_sf to match rast_com1's CRS
  nearctic_taiga_sf <- st_transform(nearctic_taiga_sf, crs_rast_com1)
}

# spatial join
joined_taiga_ne <- st_join(nearctic_taiga_sf, rast_com1, join = st_intersects)

summary(as.factor(joined_taiga_ne$sp))
#taiga, eurasian taiga, nearctic          tundra 
#1            9720            1255 
9720/nrow(joined_taiga_ne)
1255/nrow(joined_taiga_ne)
#88.6 %

#### taiga_pa ####
# Ensure nearctic_taiga is an sf object
palearctic_taiga_sf <- st_as_sf(palearctic_taiga)

if (!st_crs(palearctic_taiga_sf) == st_crs(rast_com1)) {
  # Transform nearctic_taiga_sf to match rast_com1's CRS
  palearctic_taiga_sf <- st_transform(palearctic_taiga_sf, crs_rast_com1)
}

# Perform a spatial join using st_join
joined_taiga_pa <- st_join(palearctic_taiga_sf, rast_com1, join = st_intersects)

summary(as.factor(joined_taiga_pa$sp))
#taiga, eurasian taiga, nearctic          tundra 
#13938             408            3472 
13938/nrow(joined_taiga_pa)
3472/nrow(joined_taiga_pa)
#78.2 %

#### tundra ####

# Ensure nearctic_taiga is an sf object
tundra_sf <- st_as_sf(tundra)

if (!st_crs(tundra_sf) == st_crs(rast_com1)) {
  # Transform nearctic_taiga_sf to match rast_com1's CRS
  tundra_sf <- st_transform(tundra_sf, crs_rast_com1)
}

# Perform a spatial join using st_join
joined_tundra <- st_join(tundra_sf, rast_com1, join = st_intersects)

summary(as.factor(joined_tundra$sp))
#taiga, eurasian taiga, nearctic          tundra 
#1366            1176           12337 
12337/nrow(joined_tundra)
1366/nrow(joined_tundra)
#82.9 %

(82.9 + 78.2 +88.6)/3

##### plotting ####

library(scales)
tundra <- hcl.colors(6, "Greens 3")[1] #tundra
taiga_ne <- hcl.colors(6, "Purples")[1] #nearctic
taiga_eu <- hcl.colors(6, "Blues 3")[1] #europe
col1 = c(taiga_eu, taiga_ne, tundra, "red")

plot1 <- ggplot() +
  geom_sf(data = nearctic_taiga, aes(fill = as.factor(WWF_REALM2))) +
  theme_minimal() +
  ggtitle("", subtitle = "2010")+
  labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )
print(plot1)

