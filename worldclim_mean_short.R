#worldclim mean

rm(list=ls()) #clear


#packages
library(raster)

#### tmin ####
tmin0 = list.files("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/raw/present/wc2.1_cruts4.06_2.5m_tmin_2010-2019",
             pattern = ".tif$",  
             full.names = TRUE)

tmin = stack(tmin0)

tMIN <- calc(tmin, mean) #annual average tmin

writeRaster(tMIN, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/base2015/tmin_2010_2019_mean2.grd")


