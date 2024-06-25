#1 with table 12 rows
#2 with table 4 rows, per scenario
#3 with species list and tab, 12 rows
#4 h3 with grid

#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(tidyverse)
#input

#data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays"
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

load(glue::glue("{data}/sdm_occurance_array_array_25km.rda")) #only array
load(glue::glue("{data}/sdm_occurance_array_metadata_25km.rda")) #only meta

spTable <- as.data.frame(spArrayMeta[[1]])
colnames(spTable) <- "species"
grid <- as.data.frame(spArrayMeta[[2]])

#region

ecopath <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM"

ecoreg <- st_read(glue::glue("{ecopath}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

taiga <- st_read(glue::glue("{ecopath}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga"), st_coordinates(st_centroid(.))[,2]>0)

tundra <- st_read(glue::glue("{ecopath}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Tundra"), st_coordinates(st_centroid(.))[,2]>0)

taiga_regions <- unique(taiga$ECO_NAME)
tundra_regions <- unique(tundra$ECO_NAME)

#assign tundra/taiga
grid$T <- ifelse(grid$Ecoregion %in% tundra_regions, "tundra", 
                ifelse(grid$Ecoregion %in% taiga_regions, "taiga", NA))

sp =1
species <- as.data.frame(spArray[sp,,,]) #per species
a <- as.data.frame(cbind(grid, species))
coords <- as.data.frame(st_coordinates(a$geometry, dims = c("x", "y"))) #extract coords
#Calculate decimal degree coordinates
X <- atan2(coords$Y, coords$X) * 180 / pi
Y <- atan2(coords$Y, coords$X * cos(X * pi/180)) * 180 / pi
xysp <- as.data.frame(cbind(coords, species, a$T))

colnames(xysp) <- c("x", "y", "y1s1", "y2s1", "y3s1", "y4s1", 
         "y1s2","y2s2","y3s2", "y4s2",
         "y1s3","y2s3","y3s3","y4s3", "T")

#spArray[species, cells, year, scenario]

sp1a <- as.data.frame(spArray[1,,,])

resultlist <- list()



#### base: areatab & percenttab ####

for(sp in 1:nrow(spTable)){
  
  results <- as.data.frame(matrix(ncol=12, nrow=3))
  colnames(results) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                         "y1s2","y2s2","y3s2", "y4s2",
                         "y1s3","y2s3","y3s3","y4s3")
  
  rownames(results) <- c("sum_taiga", "sum_tundra", "percent")
  
  species <- as.data.frame(spArray[sp,,,]) #per species
  a <- as.data.frame(cbind(grid, species))
  coords <- as.data.frame(st_coordinates(a$geometry, dims = c("x", "y"))) #extract coords
  xysp <- as.data.frame(cbind(coords, species, a$T))
  colnames(xysp) <- c("x", "y", "y1s1", "y2s1", "y3s1", "y4s1", 
                         "y1s2","y2s2","y3s2", "y4s2",
                         "y1s3","y2s3","y3s3","y4s3", "T")
  
  #devide
 tundra_sub <-  subset(xysp, xysp$T == "tundra")
 taiga_sub <-  subset(xysp, xysp$T == "taiga")
 
  #sum taiga
 
  for (i in 1:12) {
    v = i+2
    sum_val <- taiga_sub %>%
      pull(v) %>%  # Extract the column as a vector
      sum(na.rm = TRUE)

    results[1, i] <- sum_val
  }
 
 

 #sum tundra
 for (i in 1:12) {
   v = i+2
   sum_val <- sum(tundra_sub[, v])
   results[2, i] <- sum_val
 }

  
  
results[3,] <- results[2,]/results[1,]*100


results[1,13] <- ifelse(results[1,4] > results[1,1], "increase",
                            ifelse(results[1,4] < results[1,1], "decrease",
                                   ifelse(results[1,4] == results[1,1], "stable",
                                          NA)))

results[2,13] <- ifelse(results[2,4] > results[2,1], "increase",
                        ifelse(results[2,4] < results[2,1], "decrease",
                               ifelse(results[2,4] == results[2,1], "stable",
                                      NA)))


results[3,13] <- ifelse(results[3,4] > results[3,1], "increase",
                        ifelse(results[3,4] < results[3,1], "decrease",
                               ifelse(results[3,4] == results[3,1], "stable",
                                      NA)))


####s2

results[1,14] <- ifelse(results[1,8] > results[1,1], "increase",
                        ifelse(results[1,8] < results[1,1], "decrease",
                               ifelse(results[1,8] == results[1,1], "stable",
                                      NA)))

results[2,14] <- ifelse(results[2,8] > results[2,1], "increase",
                        ifelse(results[2,8] < results[2,1], "decrease",
                               ifelse(results[2,8] == results[2,1], "stable",
                                      NA)))


results[3,14] <- ifelse(results[3,8] > results[3,1], "increase",
                        ifelse(results[3,8] < results[3,1], "decrease",
                               ifelse(results[3,8] == results[3,1], "stable",
                                      NA)))



###s3
results[1,15] <- ifelse(results[1,12] > results[1,1], "increase",
                        ifelse(results[1,12] < results[1,1], "decrease",
                               ifelse(results[1,12] == results[1,1], "stable",
                                      NA)))

results[2,15] <- ifelse(results[2,12] > results[2,1], "increase",
                        ifelse(results[2,12] < results[2,1], "decrease",
                               ifelse(results[2,12] == results[2,1], "stable",
                                      NA)))


results[3,15] <- ifelse(results[3,12] > results[3,1], "increase",
                        ifelse(results[3,12] < results[3,1], "decrease",
                               ifelse(results[3,12] == results[3,1], "stable",
                                      NA)))


###


colnames(results)[13] <- "trend1"
colnames(results)[14] <- "trend2"
colnames(results)[15] <- "trend3"

results[3,16] <- ifelse(results[1,13] == "decrease" & results[2,13] == "decrease", "shrink",
                        ifelse(results[1,13] == "decrease" & results[2,13] == "increase", "tundra++",
                               ifelse(results[1,13] == "increase" & results[2,13] == "decrease", "taiga++",
                                      ifelse(results[1,13] == "stable" & results[2,13] == "stable", "stable",
                                             ifelse(results[1,13] == "stable" & results[2,13] == "increase", "tundra+",
                                                  ifelse(results[1,13] == "stable" & results[2,13] == "decrease", "tundra-",
                                                         ifelse(results[1,13] == "decrease" & results[2,13] == "stable", "taiga-",
                                                                ifelse(results[1,13] == "increase" & results[2,13] == "stable", "taiga+",
                                                                       ifelse(results[1,13] == "increase" & results[2,13] == "increase", "expansion", NA)))))))))

results[3,17] <- ifelse(results[1,14] == "decrease" & results[2,14] == "decrease", "shrink",
                        ifelse(results[1,14] == "decrease" & results[2,14] == "increase", "tundra++",
                               ifelse(results[1,14] == "increase" & results[2,14] == "decrease", "taiga++",
                                      ifelse(results[1,14] == "stable" & results[2,14] == "stable", "stable",
                                             ifelse(results[1,14] == "stable" & results[2,14] == "increase", "tundra+",
                                                    ifelse(results[1,14] == "stable" & results[2,14] == "decrease", "tundra-",
                                                           ifelse(results[1,14] == "decrease" & results[2,14] == "stable", "taiga-",
                                                                  ifelse(results[1,14] == "increase" & results[2,14] == "stable", "taiga+",
                                                                         ifelse(results[1,14] == "increase" & results[2,14] == "increase", "expansion", NA)))))))))

results[3,18] <- ifelse(results[1,15] == "decrease" & results[2,15] == "decrease", "shrink",
                        ifelse(results[1,15] == "decrease" & results[2,15] == "increase", "tundra++",
                               ifelse(results[1,15] == "increase" & results[2,15] == "decrease", "taiga++",
                                      ifelse(results[1,15] == "stable" & results[2,15] == "stable", "stable",
                                             ifelse(results[1,15] == "stable" & results[2,15] == "increase", "tundra+",
                                                    ifelse(results[1,15] == "stable" & results[2,15] == "decrease", "tundra-",
                                                           ifelse(results[1,15] == "decrease" & results[2,15] == "stable", "taiga-",
                                                                  ifelse(results[1,15] == "increase" & results[2,15] == "stable", "taiga+",
                                                                         ifelse(results[1,15] == "increase" & results[2,15] == "increase", "expansion", NA)))))))))


colnames(results)[16] <- "change1"
colnames(results)[17] <- "change2"
colnames(results)[18] <- "change3"

resultlist[[sp]] <- results
 
}

#### change table ####
changeT <- as.data.frame(matrix(ncol= 4, nrow=nrow(spTable)))
changeT[,1] <- spTable$species

colnames(changeT) <- c("species", "change1", "change2", "change3")
for (sp in 1:nrow(changeT)){
changeT[sp,2] <- resultlist[[sp]][3,13]
changeT[sp,3] <- resultlist[[sp]][3,14]
changeT[sp,4] <- resultlist[[sp]][3,15]
}

changeT2 <- as.data.frame(matrix(ncol= 4, nrow=nrow(spTable)))
changeT2[,1] <- spTable$species
colnames(changeT2) <- c("species", "change1", "change2", "change3")
for (sp in 1:nrow(changeT2)){
  changeT2[sp,2] <- resultlist[[sp]][3,16]
  changeT2[sp,3] <- resultlist[[sp]][3,17]
  changeT2[sp,4] <- resultlist[[sp]][3,18]
}

changeTNS <- as.data.frame(cbind(changeT$species, 
                                 changeT$change1, changeT$change2, changeT$change3,
                                 changeT2$change1, changeT2$change2, changeT2$change3,
                                 change$shift1, change$shift2, change$shift3))
colnames(changeTNS) = c("species",
                        'Taiga_Tundra1', 'Taiga_Tundra2', 'Taiga_Tundra3', 
                        "Taiga_Tundra_share1", "Taiga_Tundra_share2","Taiga_Tundra_share3",
                        "lat_minmax1", "lat_minmax2", "lat_minmax3")


changeTNS$lat_avg1 <-change_avg$summary1
changeTNS$lat_avg2 <-change_avg$summary2
changeTNS$lat_avg3 <-change_avg$summary3

setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/hypotheses_results/")
write.csv2(changeTNS, "H3_changeTNS_s3.csv", row.names = F)

summary(as.factor(changeTNS$Taiga_Tundra1))
summary(as.factor(changeTNS$Taiga_Tundra_share1))
summary(as.factor(changeTNS$lat_minmax1))
summary(as.factor(changeTNS$lat_avg1))




summary(as.factor(changeT2$change1))
summary(as.factor(changeT2$change2))
summary(as.factor(changeT2$change3))
