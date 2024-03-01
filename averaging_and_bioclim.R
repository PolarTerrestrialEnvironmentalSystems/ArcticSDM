### WorldClim data preparation for SDMs - present scenario
### cropping averaging environmental parameters (temp, precip) and calculate bioclim

library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

# wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
wd <- "/Volumes/projects/bioing/data/ArcticSDM/"

###################
#### averaging ####
###################

## Spatial extent
bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

#### files ####
fls_path <- glue::glue("{wd}raw/present")
fls      <- tibble(fls = list.files(fls_path,
                         pattern = ".tif$",  
                         recursive = TRUE)) %>%
              mutate(type  = sapply(strsplit(fls, "_"), function(x) x[4]),
                     year  = as.numeric(substr(sapply(strsplit(fls, "_"), function(x) x[length(x)]), 1, 4)),
                     month = as.numeric(substr(sapply(strsplit(fls, "_"), function(x) x[length(x)]), 6, 7)))
fls_out  <- glue::glue("{wd}environment/calibration/")


##### monthly/overall mean #####
{
  for(t in unique(fls$type)) {
    
    ### monthly
    month_list <- (fls %>% filter(type == t) %>% group_split(month)) %>% parallel::mclapply(function(m) {
      read_stars(glue::glue("{fls_path}/{m$fls}")) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
        merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames(t) %>% suppressMessages()
    }, mc.cores = 12)
    
    names(month_list) <- glue::glue("wc2.1_2.5m_{t}_2015_{1:12}")
    invisible(lapply(1:length(month_list), function(f) write_stars(month_list[[f]], glue::glue("{fls_out}monthly/{names(month_list)[f]}.tiff"))))
    
    ### annual
    overall_mean <- do.call("c", month_list) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames(t)
    write_stars(month_list[[f]], glue::glue("{fls_out}wc2.1_2.5m_{t}_2015.tiff"))
    
  }
}

##### bioclim variables #####







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
