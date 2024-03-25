#### Occurrence File Generation ####

library(megaSDM)

##### extent #####

#use extent from clipped tmin, tundra and 
#which was defined according to:!
#bbox <- st_read(glue::glue("{wd}data/Ecoregions/tnc_terr_ecoregions.shp")) %>%
#filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

tmin1 = raster("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_months1.grd")
extent(tmin1)

#extent
extent_occ <- alaska
  c(-180, 180, 44.375, 83.625)

library(rnaturalearth)

#world <- ne_countries(country = ) #save countries names
#world$sovereignt #view how names are written

usa <- ne_states(country = "United States of America") #select usa
alaska <- subset(usa, usa$iso_3166_2 == "US-AK") #select alaska

plot(alaska$geometry) #view shape
#maybe should remove little part on the other side?

tmaxalaska <- crop(tmax, alaska) #crop env data
##### species #####

tabspecies <- read.csv("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/occurrences/2023-10-23_taxa_climate_df.csv", h=T)

#use species names from ulrike/thomas-table
specieslist <- unique(tabspecies$Species)[1] #only run for 1st species now


##### download #####

occ_output <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/occurrences"

#download for these species and this extent and save on server as csv
Occurrences <- OccurrenceCollection(spplist = specieslist,
                                    output = occ_output,
                                    trainingarea = extent_occ)
