### WorldClim data preparation for SDMs - present scenario
### cropping averaging environmental parameters (temp, precip) and calculate bioclim

library(sf)
sf_use_s2(FALSE)
library(stars)
library(tibble)
library(dplyr)

# wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
# wd <- "/Volumes/projects/bioing/data/ArcticSDM/"
wd <- '/bioing/data/ArcticSDM/'

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



##############################
#### monthly/overall mean ####
##############################

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


###########################
#### bioclim variables ####
###########################

# test_bbox <- st_bbox(c(xmin = -142, xmax = -165, ymin = 54, ymax = 71), crs = 4326) %>% st_as_sfc()

fls_tif <- tibble(fls = list.files(glue::glue("{fls_out}monthly"), pattern = ".tiff")) %>%
  mutate(type  = sapply(strsplit(fls, "_"), function(x) x[[3]]),
         month = as.numeric(sapply(strsplit(fls, "_"), function(x) gsub(".tiff", "", x[[5]])))) %>%
  arrange(type, month)

{
  
  prec_tmin_tmax <- lapply(unique(fls_tif$type), function(t) {
    read_stars(glue::glue("{fls_out}monthly/{fls_tif %>% filter(type == t) %>% pull(fls)}")) %>%
      merge()
  }) %>% setNames(c("prec", "tmin", "tmax"))
  
  {        
    ## tavg
    tavg_months <- lapply(1:12, function(x) {
      (((prec_tmin_tmax$tmin[,,,x]) + (prec_tmin_tmax$tmax[,,,x])) /2) %>% adrop()
    }) %>% do.call("c", .) %>% merge()
    
    # P1. Annual Mean Temperature 
    p1 <- tavg_months %>% st_apply(., 1:2, mean) %>% setNames("p1")
    
    # P2. Mean Diurnal Range(Mean(period max-min))
    p2 <- lapply(1:12, function(x) {
      (((prec_tmin_tmax$tmin[,,,x]) - (prec_tmin_tmax$tmax[,,,x]))) %>% adrop()
    }) %>% do.call("c", .) %>% merge() %>% st_apply(., 1:2, mean) %>% setNames("p2")
    
    # P4. Temperature Seasonality (standard deviation) 
    p4 <- (tavg_months %>% st_apply(., 1:2, sd) * 100) %>% setNames("p4")
    
    # P5. Max Temperature of Warmest Period 
    p5 <- st_apply(prec_tmin_tmax$tmax, 1:2, max) %>% setNames("p5")
    
    # P6. Min Temperature of Coldest Period 
    p6 <- st_apply(prec_tmin_tmax$tmin, 1:2, min)  %>% setNames("p6")
    
    # P7. Temperature Annual Range (P5-P6) 
    p7 <- (p5 - p6) %>% setNames("p7")
    
    # P3. Isothermality (P2 / P7) 
    p3 <- ((p2 / p7) * 100) %>% setNames("p3")
    
    # P12. Annual Precipitation 
    p12 <- st_apply(prec_tmin_tmax$prec, 1:2, sum) %>% setNames("p12")
    
    # P13. Precipitation of Wettest Period 
    p13 <- st_apply(prec_tmin_tmax$prec, 1:2, max) %>% setNames("p13")
    
    # P14. Precipitation of Driest Period 
    p14 <- st_apply(prec_tmin_tmax$prec, 1:2, min) %>% setNames("p14")
    
    # P15. Precipitation Seasonality(Coefficient of Variation) 
    # the "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
    p15 <- st_apply(prec_tmin_tmax$prec+1, 1:2, raster::cv) %>% setNames("p14")
    
    # precip by quarter (3 months)		
    wet <- lapply(list(1:3, 4:6, 7:9, 10:12), function(x) {
      prec_tmin_tmax$prec[,,,x] %>% st_apply(., 1:2, sum)
    }) %>% do.call('c', .) %>% merge()
    
    # P16. Precipitation of Wettest Quarter 
    p16 <- st_apply(wet, 1:2, max) %>% setNames("p16")
    
    # P17. Precipitation of Driest Quarter 
    p17 <- st_apply(wet, 1:2, min) %>% setNames("p17")
    
    tmp <- lapply(list(1:3, 4:6, 7:9, 10:12), function(x) {
      ((tavg_months[,,,x] %>% st_apply(., 1:2, sum))/3)
    }) %>% do.call('c', .) %>% merge()
    
    # P8. Mean Temperature of Wettest Quartera
    p8 <- st_apply(c(tmp[,,,1] %>% adrop(), tmp[,,,2] %>% adrop(), tmp[,,,3] %>% adrop(), tmp[,,,4] %>% adrop(),
                     st_apply(wet, 1:2, function(x) ifelse(all(!is.na(x)), which.max(x), NA))) %>% merge(),
                   1:2, function(x) x[x[5]]) %>% setNames("p8")
    
    # P9. Mean Temperature of Driest Quarter  
    p9 <- st_apply(c(tmp[,,,1] %>% adrop(), tmp[,,,2] %>% adrop(), tmp[,,,3] %>% adrop(), tmp[,,,4] %>% adrop(),
                     st_apply(wet, 1:2, function(x) ifelse(all(!is.na(x)), which.min(x), NA))) %>% merge(),
                   1:2, function(x) x[x[5]]) %>% setNames("p9")
    
    # P10 Mean Temperature of Warmest Quarter 
    p10 <- st_apply(tmp, 1:2, max) %>% setNames("p10")
    
    # P11 Mean Temperature of Coldest Quarter
    p11 <- st_apply(tmp, 1:2, min) %>% setNames("p11")
    
    # P18. Precipitation of Warmest Quarter 
    p18 <- st_apply(c(wet[,,,1] %>% adrop(), wet[,,,2] %>% adrop(), wet[,,,3] %>% adrop(), wet[,,,4] %>% adrop(),
                      st_apply(tmp, 1:2, function(x) ifelse(all(!is.na(x)), which.max(x), NA))) %>% merge(),
                    1:2, function(x) x[x[5]]) %>% setNames("p18")
    
    # P19. Precipitation of Coldest Quarter 
    p19 <- st_apply(c(wet[,,,1] %>% adrop(), wet[,,,2] %>% adrop(), wet[,,,3] %>% adrop(), wet[,,,4] %>% adrop(),
                      st_apply(tmp, 1:2, function(x) ifelse(all(!is.na(x)), which.min(x), NA))) %>% merge(),
                    1:2, function(x) x[x[5]]) %>% setNames("p19")
    
    
    biovar <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                p11, p12, p13, p14, p15, p16, p17, p18, p19) %>% merge() %>% setNames("biovar")
    
    
    write_stars(biovar, glue::glue("{fls_out}wc2.1_2.5m_biovar_2015.tiff"))
  }
  
  }