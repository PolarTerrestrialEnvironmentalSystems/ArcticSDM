#monthly tmin

library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
#wd <- "/Volumes/projects/bioing/data/ArcticSDM/"

## We only need tundra and taiga
bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)


### tmin
fls_tmin <- list.files(glue::glue("{wd}raw/present/wc2.1_cruts4.06_2.5m_tmin_2010-2019"),
                       pattern = ".tif$",  
                       full.names = TRUE)

fls_month  <- as.numeric(substr(sapply(strsplit(fls_tmin, "_"), function(x) x[length(x)]), 6, 7))

## first do it by month
tmin_stars_months <- lapply(unique(fls_month), function(x) {
  read_stars(fls_tmin[fls_month==x]) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
    merge() %>% st_apply(., 1:2, mean, future = T) %>% suppressMessages()
})
a <- unlist(tmin_stars_months)
stack(tmin_stars_months)
library(raster)
raster(tmin_stars_months)
writeRaster(tmin_stars_months, 'tmin_stars_months.tif')
print(tmin_stars_months)
sink()
b <- tmin_stars_months[[1]]
write_stars(tmin_stars_months[[1]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months1.grd")
write_stars(tmin_stars_months[[2]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months2.grd")
write_stars(tmin_stars_months[[3]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months3.grd")
write_stars(tmin_stars_months[[4]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months4.grd")
write_stars(tmin_stars_months[[5]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months5.grd")
write_stars(tmin_stars_months[[6]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months6.grd")
write_stars(tmin_stars_months[[7]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months7.grd")
write_stars(tmin_stars_months[[8]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months8.grd")
write_stars(tmin_stars_months[[9]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months9.grd")
write_stars(tmin_stars_months[[10]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months10.grd")
write_stars(tmin_stars_months[[11]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months11.grd")
write_stars(tmin_stars_months[[12]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months12.grd")

write_stars(tmin_stars_months[[1]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months1")
write_stars(tmin_stars_months[[2]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months2")
write_stars(tmin_stars_months[[3]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months3")
write_stars(tmin_stars_months[[4]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months4")
write_stars(tmin_stars_months[[5]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months5")
write_stars(tmin_stars_months[[6]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months6")
write_stars(tmin_stars_months[[7]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months7")
write_stars(tmin_stars_months[[8]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months8")
write_stars(tmin_stars_months[[9]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months9")
write_stars(tmin_stars_months[[10]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months10")
write_stars(tmin_stars_months[[11]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months11")
write_stars(tmin_stars_months[[12]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months12")

tmin <- do.call("c", tmin_stars_months) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames("tmin_2015")
plot(tmin)

write_stars(tmin, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_mean3")
