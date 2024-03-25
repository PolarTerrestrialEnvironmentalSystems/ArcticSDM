#monthly



library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
#wd <- "/Volumes/projects/bioing/data/ArcticSDM/"

## We only need tundra and taiga
bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)


### prec
fls_prec <- list.files(glue::glue("{wd}raw/present/wc2.1_cruts4.06_2.5m_prec_2010-2019"),
                       pattern = ".tif$",  
                       full.names = TRUE)

fls_month  <- as.numeric(substr(sapply(strsplit(fls_prec, "_"), function(x) x[length(x)]), 6, 7))

## first do it by month
prec_stars_months <- lapply(unique(fls_month), function(x) {
  read_stars(fls_prec[fls_month==x]) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
    merge() %>% st_apply(., 1:2, mean, future = T) %>% suppressMessages()
})

prec <- do.call("c", prec_stars_months) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames("prec_2015")
plot(prec)

write_stars(prec, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_mean")
