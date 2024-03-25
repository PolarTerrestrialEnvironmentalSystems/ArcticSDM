#compare bioclims with dismo and by Simeon

library(terra)
bio_simeon = rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/wc2.1_2.5m_biovar_2015.tiff")
plot(bio_simeon$p1)
bio_simeon

bio_dismo = rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/calibration/wc2.1_2.5m_biovar_2015.tiff")
plot(bio_dismo$p1)
bio_dismo