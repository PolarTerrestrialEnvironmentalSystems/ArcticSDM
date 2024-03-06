# install.packages("rgbif")
library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

## Spatial extent
ecoreg <- st_read(glue::glue("/Users/slisovsk/Google Drive/My Drive/GeoDat/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

tmp_folder <- dir.create("")
save_dir   <- ""

for(i in 1:nrow(ecoreg)) {
  
  if(!file.exists(glue::glue("/Users/slisovsk/Desktop/{ecoreg[i,] %>% pull('ECO_ID_U')}_gbif_all.csv"))) {
    large_wkt <- ecoreg[i,] %>% 
      st_bbox() %>% st_as_sfc(crs = 4326) %>%
      sf::st_as_text()
    
    file_list <- occ_download(pred_within(large_wkt),format = "SIMPLE_CSV")
    occ_download_wait(file_list[1])
    
    dwnl <- occ_download_get(file_list[1], overwrite = T) %>%
      occ_download_import(path = tmp_folder)
    
    inPoly <- dwnl %>% st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
      st_intersection(., ecoreg[i,] %>% st_geometry()) %>% 
      suppressMessages() %>% suppressWarnings()

    write_csv(inPoly %>% st_drop_geometry(), 
              file = glue::glue("/Users/slisovsk/Desktop/{ecoreg[i,] %>% pull('ECO_ID_U')}_gbif_all.csv"))
  }

}

