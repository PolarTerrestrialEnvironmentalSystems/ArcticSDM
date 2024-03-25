#worldclim mean

rm(list=ls()) #clear


#packages
library(raster)

#### tmin ####
tmax0 = list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/raw/present/wc2.1_cruts4.06_2.5m_tmax_2010-2019",
             pattern = ".tif$",  
             full.names = TRUE)

tmax = stack(tmax0)

tMAX <- calc(tmax, mean) #annual average tmin

writeRaster(tMAX, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmax_2010_2019_mean.grd")

