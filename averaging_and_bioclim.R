#averaging parameters and calculate bioclim
#still to redundant, needs to be coded smarter

#### averaging ####

#preparations
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"

## We only need tundra and taiga
bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

#### tmin ####

{
  fls_tmin <- list.files(glue::glue("{wd}raw/present/wc2.1_cruts4.06_2.5m_tmin_2010-2019"),
                       pattern = ".tif$",  
                       full.names = TRUE)

fls_month  <- as.numeric(substr(sapply(strsplit(fls_tmin, "_"), function(x) x[length(x)]), 6, 7))

##### monthly #####

# calculate monthly
tmin_stars_months <- lapply(unique(fls_month), function(x) {
  read_stars(fls_tmin[fls_month==x]) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
    merge() %>% st_apply(., 1:2, mean, future = T) %>% suppressMessages()
})

# export monthly
write_stars(tmin_stars_months[[1]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months1.grd")
write_stars(tmin_stars_months[[2]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months2.grd")
write_stars(tmin_stars_months[[3]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months3.grd")
write_stars(tmin_stars_months[[4]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months4.grd")
write_stars(tmin_stars_months[[5]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months5.grd")
write_stars(tmin_stars_months[[6]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months6.grd")
write_stars(tmin_stars_months[[7]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months7.grd")
write_stars(tmin_stars_months[[8]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months8.grd")
write_stars(tmin_stars_months[[9]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months9.grd")
write_stars(tmin_stars_months[[10]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months10.grd")
write_stars(tmin_stars_months[[11]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months11.grd")
write_stars(tmin_stars_months[[12]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_months12.grd")

##### yearly #####
tmin <- do.call("c", tmin_stars_months) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames("tmin_2015")

write_stars(tmin, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmin_2010_2019_mean3")

}

#### tmax ####

{
fls_tmax <- list.files(glue::glue("{wd}raw/present/wc2.1_cruts4.06_2.5m_tmax_2010-2019"),
                       pattern = ".tif$",  
                       full.names = TRUE)

fls_month  <- as.numeric(substr(sapply(strsplit(fls_tmax, "_"), function(x) x[length(x)]), 6, 7))

##### monthly #####

# calculate monthly
tmax_stars_months <- lapply(unique(fls_month), function(x) {
  read_stars(fls_tmax[fls_month==x]) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
    merge() %>% st_apply(., 1:2, mean, future = T) %>% suppressMessages()
})

# export monthly
write_stars(tmax_stars_months[[1]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months1.grd")
write_stars(tmax_stars_months[[2]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months2.grd")
write_stars(tmax_stars_months[[3]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months3.grd")
write_stars(tmax_stars_months[[4]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months4.grd")
write_stars(tmax_stars_months[[5]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months5.grd")
write_stars(tmax_stars_months[[6]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months6.grd")
write_stars(tmax_stars_months[[7]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months7.grd")
write_stars(tmax_stars_months[[8]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months8.grd")
write_stars(tmax_stars_months[[9]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months9.grd")
write_stars(tmax_stars_months[[10]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months10.grd")
write_stars(tmax_stars_months[[11]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months11.grd")
write_stars(tmax_stars_months[[12]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_months12.grd")

##### yearly #####
tmax <- do.call("c", tmax_stars_months) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames("tmax_2015")

write_stars(tmax, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/tmax_2010_2019_mean3")

}

#### prec ####

{
  
fls_prec <- list.files(glue::glue("{wd}raw/present/wc2.1_cruts4.06_2.5m_prec_2010-2019"),
                       pattern = ".tif$",  
                       full.names = TRUE)

fls_month  <- as.numeric(substr(sapply(strsplit(fls_prec, "_"), function(x) x[length(x)]), 6, 7))

##### monthly #####

# calculate monthly
prec_stars_months <- lapply(unique(fls_month), function(x) {
  read_stars(fls_prec[fls_month==x]) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
    merge() %>% st_apply(., 1:2, mean, future = T) %>% suppressMessages()
})

# export monthly
write_stars(prec_stars_months[[1]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months1.grd")
write_stars(prec_stars_months[[2]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months2.grd")
write_stars(prec_stars_months[[3]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months3.grd")
write_stars(prec_stars_months[[4]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months4.grd")
write_stars(prec_stars_months[[5]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months5.grd")
write_stars(prec_stars_months[[6]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months6.grd")
write_stars(prec_stars_months[[7]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months7.grd")
write_stars(prec_stars_months[[8]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months8.grd")
write_stars(prec_stars_months[[9]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months9.grd")
write_stars(prec_stars_months[[10]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months10.grd")
write_stars(prec_stars_months[[11]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months11.grd")
write_stars(prec_stars_months[[12]], "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_months12.grd")

##### yearly #####
prec <- do.call("c", prec_stars_months) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames("prec_2015")

write_stars(prec, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/prec_2010_2019_mean3")

}

#### bioclim ####

{
library(raster)

##### prepare tmin #####
tmin1 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months1.grd")
tmin2 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months2.grd")
tmin3 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months3.grd")
tmin4 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months4.grd")
tmin5 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months5.grd")
tmin6 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months6.grd")
tmin7 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months7.grd")
tmin8 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months8.grd")
tmin9 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months9.grd")
tmin10 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months10.grd")
tmin11 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months11.grd")
tmin12 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months12.grd")


tmin2015 = stack(tmin1, tmin2, tmin3, tmin4, tmin5, tmin6, tmin7, tmin8, tmin9, 
                  tmin10, tmin11, tmin12)

month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

names(tmin2015) <- month

##### prepare tmax #####
tmax1 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months1.grd")
tmax2 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months2.grd")
tmax3 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months3.grd")
tmax4 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months4.grd")
tmax5 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months5.grd")
tmax6 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months6.grd")
tmax7 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months7.grd")
tmax8 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months8.grd")
tmax9 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months9.grd")
tmax10 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months10.grd")
tmax11 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months11.grd")
tmax12 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_months12.grd")


tmax2015 = stack(tmax1, tmax2, tmax3, tmax4, tmax5, tmax6, tmax7, tmax8, tmax9, 
                 tmax10, tmax11, tmax12)

month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

names(tmax2015) <- month

##### prepare prec #####
prec1 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months1.grd")
prec2 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months2.grd")
prec3 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months3.grd")
prec4 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months4.grd")
prec5 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months5.grd")
prec6 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months6.grd")
prec7 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months7.grd")
prec8 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months8.grd")
prec9 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months9.grd")
prec10 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months10.grd")
prec11 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months11.grd")
prec12 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/prec_2010_2019_months12.grd")


prec2015 = stack(prec1, prec2, prec3, prec4, prec5, prec6, prec7, prec8, prec9, 
                 prec10, prec11, prec12)

month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

names(prec2015) <- month

##### calculate #####
library(dismo)

bio2015 = biovars(prec = prec2015, tmin = tmin2015, tmax = tmax2015)

writeRaster(bio2015, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/server2015/bio2015.grd")

}
