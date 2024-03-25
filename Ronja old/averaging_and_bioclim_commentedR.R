### WorldClim data preparation for SDMs - present scenario
### cropping averaging environmental parameters (temp, precip) and calculate bioclim

library(sf)
sf_use_s2(FALSE)
library(stars)
library(tibble)
library(dplyr)

wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/"
# wd <- "/Volumes/projects/bioing/data/ArcticSDM/"
#wd <- '/bioing/data/ArcticSDM/'

###################
#### averaging ####
###################

## Spatial extent
bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)


#glue: Expressions enclosed by braces will be evaluated as R code, here: wd
#st_coordinates: retrieve coordinates for raster or vector cube cells
  #matrix with coordinates (X, Y, possibly Z and/or M) in rows,
  #possibly followed by integer indicators L1,...,L3 that point out to which structure the coordinate belongs; 
  #for POINT this is absent (each coordinate is a feature), 
#st_centroid() is used to extract the centroid coordinates of each polygon feature and returns a point geometry. 
  #Note that the centroid is not always inside of the polygon that it is the center of 
  #(for example, a “C-shaped” island or a doughnut).
#st_coordinates(st_centroid(.))[,2]>0) wollen Y = latidude > 0 filtern

#### files ####
fls_path <- glue::glue("{wd}raw/present")
fls      <- tibble(fls = list.files(fls_path,
                                    pattern = ".tif$",  
                                    recursive = TRUE)) %>%
  mutate(type  = sapply(strsplit(fls, "_"), function(x) x[4]),
         year  = as.numeric(substr(sapply(strsplit(fls, "_"), function(x) x[length(x)]), 1, 4)),
         month = as.numeric(substr(sapply(strsplit(fls, "_"), function(x) x[length(x)]), 6, 7)))
fls_out  <- glue::glue("{wd}environment/calibration/")

#recursive: logical. Should the listing recurse into directories?
#function(x): an anonymous function. The function basically returns the same objects (= does nothing). 
  #x[4] returns the 4th element of the splitted elements
  #x[length(x)]) returns the last element of the splitted elements
#substr(x, start, stop) 1,4 from 1 to 4 e.g.  digits 2015 


##############################
#### monthly/overall mean ####
##############################

{
  for(t in unique(fls$type)) { 
    #for prec, tmax, tmin
    #filter for these (counts with t not i)
    
    ### monthly
    #calculate
    month_list <- (fls %>% filter(type == t) %>% group_split(month)) %>% parallel::mclapply(function(m) {
      read_stars(glue::glue("{fls_path}/{m$fls}")) %>% st_crop(bbox %>% st_bbox() %>% st_as_sfc(crs = 4326)) %>% 
        merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames(t) %>% suppressMessages()
    }, mc.cores = 12)
    
    #parallel::mclapply is a parallelized version of lapply, it returns a list of the same length as X, 
      #each element of which is the result of applying FUN to the corresponding element of X. 
      #It relies on forking and hence is not available on Windows unless mc.cores = 1. 
    #st_apply apply a function to array dimensions: aggregate over space, time, or something else
    #MARGIN see apply; index number(s) or name(s) of the dimensions over which FUN will be applied
    #FUTURE 	logical;if TRUE, use future.apply::future_apply
    
    #write
    names(month_list) <- glue::glue("wc2.1_2.5m_{t}_2015_{1:12}")
    invisible(lapply(1:length(month_list), function(f) write_stars(month_list[[f]], 
                      glue::glue("{fls_out}monthly/{names(month_list)[f]}.tiff"))))
    
    #invisible: Return a (temporarily) invisible copy of an object.
      #This function can be useful when it is desired to have functions return values which can be assigned, 
      #but which do not print when they are not assigned. 
    #f loops from 1:12 through the months
    #glue glues the name and path of the saved files together
    
    ### annual
    overall_mean <- do.call("c", month_list) %>% merge() %>% st_apply(., 1:2, mean, future = T) %>% setNames(t)
    write_stars(month_list[[f]], glue::glue("{fls_out}wc2.1_2.5m_{t}_2015.tiff"))
    
    #do.call: constructs and executes a function call from a name or a function and a list of arguments to be passed to it.
    #t is here again the variable prec tmax tmin
  }
}


###########################
#### bioclim variables ####
###########################

# test_bbox <- st_bbox(c(xmin = -142, xmax = -165, ymin = 54, ymax = 71), crs = 4326) %>% st_as_sfc()
#a numeric vector of length four, with xmin, ymin, xmaxand ymax values; 
#if obj is of class sf, sfc, Spatial or Raster, the object returned has a class bbox, 
#an attribute crs and a method to print the bbox 
#and an st_crs method to retrieve the coordinate reference system corresponding to obj d
#(and hence the bounding box). 
#st_as_sfc has a methods for bbox objects to generate a polygon around the four bounding box points.

fls_tif <- tibble(fls = list.files(glue::glue("{fls_out}monthly"), pattern = ".tiff")) %>%
  mutate(type  = sapply(strsplit(fls, "_"), function(x) x[[3]]),
         month = as.numeric(sapply(strsplit(fls, "_"), function(x) gsub(".tiff", "", x[[5]])))) %>%
  arrange(type, month)

#x[[3]] prec/tmax/tmin
#x[[5]] month

{
  
  #sorts files by variable
  prec_tmin_tmax <- lapply(unique(fls_tif$type), function(t) {
    read_stars(glue::glue("{fls_out}monthly/{fls_tif %>% filter(type == t) %>% pull(fls)}")) %>%
      merge()
  }) %>% setNames(c("prec", "tmin", "tmax"))
  
  #pull:is similar to $, selects a column in a data frame and transforms it into a vector
  
  {        
    ## tavg
    tavg_months <- lapply(1:12, function(x) {
      (((prec_tmin_tmax$tmin[,,,x]) + (prec_tmin_tmax$tmax[,,,x])) /2) %>% adrop()
    }) %>% do.call("c", .) %>% merge()
    
    #prec_tmin_tmax$tmin[,,,x] loop durch x = month, Rest wird genommen
    #atavg = tmin + tmax/ 2
    #adrop: Drop degenerate dimensions of an array object. 
    
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
    p15 <- st_apply(prec_tmin_tmax$prec+1, 1:2, raster::cv) %>% setNames("p15")
    
    #raster::cv Compute the coefficient of variation (expressed as a percentage).
    
    # precip by quarter (3 months)		
    wet <- lapply(list(1:3, 4:6, 7:9, 10:12), function(x) {
      prec_tmin_tmax$prec[,,,x] %>% st_apply(., 1:2, sum)
    }) %>% do.call('c', .) %>% merge()
    
    # P16. Precipitation of Wettest Quarter 
    p16 <- st_apply(wet, 1:2, max) %>% setNames("p16")
    
    # P17. Precipitation of Driest Quarter 
    p17 <- st_apply(wet, 1:2, min) %>% setNames("p17")
    
    #mean temperature of quarters
    tmp <- lapply(list(1:3, 4:6, 7:9, 10:12), function(x) {
      ((tavg_months[,,,x] %>% st_apply(., 1:2, sum))/3)
    }) %>% do.call('c', .) %>% merge()
    
    #/3: divided by three months of quarter
    
    # P8. Mean Temperature of Wettest Quartera
    p8 <- st_apply(c(tmp[,,,1] %>% adrop(), tmp[,,,2] %>% adrop(), tmp[,,,3] %>% adrop(), tmp[,,,4] %>% adrop(),
                     st_apply(wet, 1:2, function(x) ifelse(all(!is.na(x)), which.max(x), NA))) %>% merge(),
                   1:2, function(x) x[x[5]]) %>% setNames("p8")
    #tmp[,,,1] %>% adrop()... call all quarters and remove dimensions
    #which.max returns the position of the element with the maximal value in a vector. 
    #x[x[5]] indicates the max value which was identified and added as column
    
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