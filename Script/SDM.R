library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)


# wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
wd <- "/Volumes/projects/bioing/data/ArcticSDM/"
# wd <- '/bioing/data/ArcticSDM/'


### Grid
grid <- st_read(glue::glue("{wd}/Data/occurance_grid.shp"))

### Shoreline map
map <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf') %>%
  st_intersection(grid %>% st_bbox(crs = 4326) %>% st_as_sfc())


