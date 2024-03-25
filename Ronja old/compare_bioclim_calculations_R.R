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
bio_rcp452030 = rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/RCP45/2030/wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040.tif")
bio_rcp = crop(bio_rcp452030, bio_simeon)
bio_simeon
bio_rcp 

par(mfrow = c(1, 2)) 
plot(bio_simeon$p1)
plot(bio_rcp$`wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_1`)
hist(bio_simeon$p3)
hist(bio_rcp$`wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_12`)
summary(bio_simeon$p3)
summary(bio_rcp$`wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_3`)
