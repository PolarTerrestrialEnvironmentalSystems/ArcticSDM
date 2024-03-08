#bioclim

rm(list=ls()) #clear
getwd()

library(raster)

month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

#tmin2020
tmin2020 = list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/monthly/tmin",
                      pattern = ".tiff$",  
                      full.names = TRUE)


tmin = stack(tmin2020)

names(tmin) <- month

#tmax2020
tmax2020 = list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/monthly/tmax",
                      pattern = ".tiff$",  
                      full.names = TRUE)

tmax =stack(tmax2020)

names(tmax) <- month


#prec2020
prec2020 = list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/monthly/prec",
                      pattern = ".tiff$",  
                      full.names = TRUE)

prec =stack(prec2020)

names(prec) <- month

#biovars

library(dismo)

bio2020 = biovars(prec = prec, tmin = tmin, tmax = tmax)

writeRaster(bio2020, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/wc2.1_2.5m_biovar_2015_dismo.grd")


#compare bioclims with dismo and by Simeon

bio_dismo = bio2020
#bio_dismo = rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/wc2.1_2.5m_biovar_2015_dismo.grd")
plot(bio_dismo$p1)
bio_dismo

library(terra)
bio_simeon = rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/wc2.1_2.5m_biovar_2015.tiff")
plot(bio_simeon$p1)
bio_simeon

