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

path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

### Maps
          
ecoreg   <- st_read(glue::glue("{path}/tnc_terr_ecoregions.shp")) %>%
              filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

load(glue::glue("{path}/occTab_all.rda"))

## reduce occTab
minSpecies <- 150
occTab <- occTab_all %>% filter(species %in% (occTab_all %>% group_by(species) %>% st_drop_geometry() %>% summarise(count = n()) %>%
                                                arrange(desc(count)) %>% filter(count>minSpecies) %>% pull(species))) %>%
  mutate(ecoreg = (ecoreg %>% pull(ECO_NAME))[apply(st_intersects(., ecoreg %>% dplyr::select(ECO_NAME) %>% st_transform(proj), sparse = FALSE), 1, function(x) which(x)[1])],
         .after = 'species')


# Define the polar projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"


# Get world map data and transform it to WGS84
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, crs = 4326)

#### occplot ####
occTab$sp = 1
occTab$sp <- as.factor(occTab$sp)

occs <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "white") +  # Background map
  geom_sf(data = occTab, size = 0.1, aes(color = sp),show.legend = FALSE) +  # Keep mapping to sp
  scale_color_manual(values = rep("black", length(unique(occTab$sp)))) +  # Set all species colors to black
  coord_sf(crs = proj, 
           xlim = c(-5000000, 5000000),  # Adjust these values based on your projection
           ylim = c(-5000000, 4000000),  # Adjust these values based on your projection
           expand = FALSE) +
  theme_minimal()+
  labs(title = "A)
       ")
  

print(occs)

#### with p1B_sankey script ####
library(cowplot)
combined <- plot_grid(occs, sankey_plot, ncol = 2)
combined <- plot_grid(occs, sankey_plot, ncol = 2, 
                      #labels = c("A)", "B)"), 
                      label_size = 14)
print(combined)


