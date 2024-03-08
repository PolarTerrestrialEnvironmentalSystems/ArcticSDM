#### packages ####
library(SSDM)
library(raster)

#### input ####

#test data input and SSDM can be activated, is commented out

##### env #####

files = list.files("worldclim/means",
                   pattern = ".grd$",  
                   full.names = TRUE)

Env = stack(files)


##### occurrences #####

#Data from gbif-table
# Occ = load_occ(path = getwd(), Env, file = 'data/speciesgbif3000.csv',
#                Xcol = 'x', Ycol = 'y', Spcol = 'species', GeoRes = TRUE, sep = ',',
#                verbose = FALSE, GUI = FALSE)


Picea_glauca = read.csv("occurrences3/Picea_glauca.csv")
Achillea_millefolium = read.csv("occurrences3/Achillea_millefolium.csv")
Chamaenerion_angustifolium = read.csv("occurrences3/Chamaenerion_angustifolium.csv")

speciesdata = rbind(Picea_glauca, Achillea_millefolium, Chamaenerion_angustifolium)

#### SSDM ####

SSDM_GAM3 <- stackmod(c("GAM"), speciesdata, 
                      Env, rep = 1, Xcol = 'x', Ycol = 'y',
                      ensemble.metric = c('AUC'), weight = TRUE,
                      ensemble.thresh = 0, Pcol=NULL,tmp=FALSE,
                      Spcol = 'Species', verbose = FALSE)