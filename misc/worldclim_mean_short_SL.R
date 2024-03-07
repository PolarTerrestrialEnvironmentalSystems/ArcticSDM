library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

# wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
wd <- "/Volumes/projects/bioing/data/ArcticSDM/"

## We only need tundra and taiga
bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)


### tmin
fls_tmin <- list.files(glue::glue("{wd}raw/present/wc2.1_cruts4.06_2.5m_tmin_2010-2019"),
                        pattern = ".tif$",  
                        full.names = TRUE)

fls_year  <- as.numeric(substr(sapply(strsplit(fls_tmin, "_"), function(x) x[length(x)]), 1, 4))

## first do it by year or month
tmin_stars_years <- lapply(unique(fls_year), function(x) {
  read_stars(fls_tmin[fls_year==x]) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
    merge() %>% st_apply(., 1:2, mean, future = T) %>% suppressMessages()
})

tmin <- do.call("c", tmin_stars_years) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames("tmin_2015")
plot(tmin)

write_stars(tmin, "...")


