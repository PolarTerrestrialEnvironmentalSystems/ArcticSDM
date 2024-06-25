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

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays"

load(glue::glue("{data}/spTable.rda"))
load(glue::glue("{data}/binaryArray_25km.rda"))
load(glue::glue("{data}/grid_25km.rda"))
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
species <- as.data.frame(binarArray[sp,,,]) #per species
a <- as.data.frame(cbind(grid, species))
coords <- as.data.frame(st_coordinates(a$geometry, dims = c("x", "y"))) #extract coords
xysp <- as.data.frame(cbind(coords, species, a$T))

colnames(xysp) <- c("x", "y", "y1s1", "y2s1", "y3s1", "y4s1", 
         "y1s2","y2s2","y3s2", "y4s2",
         "y1s3","y2s3","y3s3","y4s3", "T")

#binarArray[species, cells, year, scenario]

sp1a <- as.data.frame(binarArray[1,,,])

resultlist <- list()

results <- as.data.frame(matrix(ncol=14, nrow=3))
colnames(results) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                       "y1s2","y2s2","y3s2", "y4s2",
                       "y1s3","y2s3","y3s3","y4s3", "trend", "change")

rownames(results) <- c("sum_taiga", "sum_tundra", "percent")

#### base: areatab & percenttab ####

for(sp in 1:nrow(spTable)){
  
  sp = 1
  species <- as.data.frame(binarArray[sp,,,]) #per species
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
  for (i in 1:ncol(results)) {
    v = i+2
    sum_val <- sum(taiga_sub[, v])
    taiga_sub[1, 3]
    results[1, i] <- sum_val
  }
 
 

 #sum tundra
 for (i in 1:ncol(results)) {
   v = i+2
   sum_val <- sum(tundra_sub[, v])
   results[2, i] <- sum_val
 }
}
  
  
results[3,] <- results[2,]/results[1,]*100


results[1,13] <- ifelse(mean(results[1,4], results[1,8], results[1,12]) > results[1,1], "increase",
                            ifelse(mean(results[1,4], results[1,8], results[1,12]) < results[1,1], "decrease",
                                   ifelse(mean(results[1,4], results[1,8], results[1,12]) == results[1,1], "stable",
                                          NA)))

results[2,13] <- ifelse(mean(results[2,4], results[2,8], results[2,12]) > results[2,1], "increase",
                        ifelse(mean(results[2,4], results[2,8], results[2,12]) < results[2,1], "decrease",
                               ifelse(mean(results[2,4], results[2,8], results[2,12]) == results[2,1], "stable",
                                      NA)))

results[3,13] <- ifelse(mean(results[3,4], results[3,8], results[3,12]) > results[3,1], "tundra+",
                        ifelse(mean(results[3,4], results[3,8], results[3,12]) < results[3,1], "taiga+",
                               ifelse(mean(results[3,4], results[3,8], results[3,12]) == results[3,1], "stable",
                                      NA)))

results[3,14] <- ifelse(results[1,13] == "decrease" & results[2,13] == "decrease", "shrink",
                        ifelse(results[1,13] == "decrease" & results[2,13] == "increase", "tundra+",
                               ifelse(results[1,13] == "increase" & results[2,13] == "decrease", "taiga+",
                                      ifelse(results[1,13] == "stable" & results[2,13] == "stable", "stable",
                                             ifelse(results[1,13] == "increase" & results[2,13] == "increase", "expansion", NA)))))













#########################################
sp1a <- as.data.frame(binarArray[1,,,]) #per species

sp <- as.data.frame((sp1a))
#tsp <- as.data.frame(t(sp1a))

a <- as.data.frame(cbind(grid, sp))
#a <- as.data.frame(cbind(grid, tsp))
stg

coords <- as.data.frame(st_coordinates(a$geometry, dims = c("x", "y"))) #extract coords

b <- as.data.frame(cbind(coords, sp))

stg <- st_as_stars(grid$geometry, breaks="equal")
stg <- st_as_stars(grid)

plot(stg)
stb
########################################
stg <- st_as_stars(grid)
tgrid <- as.data.frame(t(grid))

gridsp <- c(grid$geometry, sp1a)

g <- stg %>% 
  mutate(value = sp)

g <- st_make_grid(grid$geometry) %>% st_centroid()

tgrid[5,] <- NA
sp1a[1,]

plot(grid$geometry)

grid2 <- gridsp %>%
  st_coordinates(geometry) %>%
  st_as_sf(.) 


plot(gridsp)

gridsp

resultlist <- list()
percenttab <- as.data.frame(matrix(nrow = nrow(spTable), ncol=14))
colnames(percenttab) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                       "y1s2","y2s2","y3s2", "y4s2",
                       "y1s3","y2s3","y3s3","y4s3", "trend", "species")

areatab <- as.data.frame(matrix(nrow = nrow(spTable), ncol=13))
colnames(areatab) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                          "y1s2","y2s2","y3s2", "y4s2",
                          "y1s3","y2s3","y3s3","y4s3", "species")

#### base: areatab & percenttab ####
for(sp in 1:nrow(spTable)){


sp1a <- as.data.frame(binarArray[sp,,,])

results <- as.data.frame(matrix(ncol=12, nrow=3))

colnames(results) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                       "y1s2","y2s2","y3s2", "y4s2",
                       "y1s3","y2s3","y3s3","y4s3")

rownames(results) <- c("sum", "delta", "percent")


#sum
for (i in 1:ncol(results)) {
  # Calculate the sum of values in the first 38657 rows for the current column
  sum_val <- sum(sp1a[1:38657, i])
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[1, i] <- sum_val
}



#delta
for (i in 1:3) {
  # Calculate the sum of values in the first 38657 rows for the current column
  delta_val <- results[1, i+1] - results[1, 1]
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[2, i+1] <- delta_val
}

for (i in 1:3) {
  # Calculate the sum of values in the first 38657 rows for the current column
  delta_val <- results[1, i+5] - results[1, 1]
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[2, i+5] <- delta_val
}

for (i in 1:3) {
  # Calculate the sum of values in the first 38657 rows for the current column
  delta_val <- results[1, i+9] - results[1, 1]
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[2, i+9] <- delta_val
}

#percent
for (i in 1:ncol(results)) {
  # Calculate the sum of values in the first 38657 rows for the current column
  per_val <- results[2,i]/results[1,1]*100
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[3, i] <- per_val
}

#resultlist
resultlist[[sp]] <- results

#### percenttab H0 ####
percenttab[sp,] <- results[3,]
percenttab[sp,14] <- spTable[sp, 1]

percenttab[sp,13] <- ifelse(mean(results[3,4], results[3,8], results[3,12]) > 0, "increase",
                        ifelse(mean(results[3,4], results[3,8], results[3,12]) < 0, "decrease",
                               ifelse(mean(results[3,4], results[3,8], results[3,12]) == 0, "stable",
                                      NA)))
percenttab[sp, 15] <- ifelse(results[3,4]> 0, "increase",
                             ifelse(results[3,4]< 0, "decrease",
                                    ifelse(results[3,4] == 0, "stable",
                                           NA)))
percenttab[sp, 16] <- ifelse(results[3,8]> 0, "increase",
                             ifelse(results[3,8]< 0, "decrease",
                                    ifelse(results[3,8] == 0, "stable",
                                           NA)))
percenttab[sp, 17] <- ifelse(results[3,12]> 0, "increase",
                             ifelse(results[3,12]< 0, "decrease",
                                    ifelse(results[3,12] == 0, "stable",
                                           NA)))
colnames(percenttab) [15:17] <- c("trend_s1", "trend_s2", "trend_s3")
#areatab H1
areatab[sp,] <- results[1,]
areatab[sp,13] <- spTable[sp, 1]

}


summary(as.factor(percenttab$trend))

#### gform H1 based on areatab ####
traits <- read.csv("C:/Users/roschw001/Documents/migration values/modspeciesdist.csv", sep=",", dec=",")
traits_gform <- as.data.frame(cbind(traits$species, traits$gform))
colnames(traits_gform) <- c("species", "gform")
areatab <- areatab %>% left_join(traits_gform, by="species") 

summary(as.factor(gform_result$gform))
gform_result <- as.data.frame(matrix(nrow=7, ncol=14))  
colnames(gform_result) = colnames(areatab)
gform_result$gform <-  unique(traits_gform$gform)

##### gformarea #####
library(dplyr)

#sum area for all yxsx
r <- seq(1,13,1)
areatab <- rbind(r, areatab)

areatab[1,13] <- "sum"
for (i in 1:12){
  # Calculate the sum of valuesin the first 38657 rows for the current column
  sum_val <- sum(areatab[2:17, i])
  areatab[1, i] <- sum_val
}

gformarea <- areatab %>%
  group_by(gform) %>%
  summarise(across(1:12, ~ sum(.x, na.rm = TRUE)))
areatab %>% filter(gform)


for (i in 2:13){
gformarea[6,i] <- sum(gformarea[1:5, i])
}
gformarea[6,1] = 'all'

##### gformpercent #####
gformpercent = gformarea

for (n in 1:nrow(gformpercent)) {
  for (i in 2:13) {
    
    
    gformpercent[n,i] <- gformarea[n, i]/gformarea[6, i]*100
  }}

##### trend #####
#trend s1
for (i in 1:nrow(gformpercent)){
  gformpercent[i,14] <- ifelse(gformpercent[i,2] < gformpercent[i,5], "increase2100",
                         ifelse(gformpercent[i,2] > gformpercent[i,5], "decrease2100",
                         ifelse(gformpercent[i,2] == gformpercent[i,5], "stable", NA)))
          
}

#trend s2
for (i in 1:nrow(gformpercent)){
  gformpercent[i,15] <- ifelse(gformpercent[i,6] < gformpercent[i,9], "increase2100",
                               ifelse(gformpercent[i,6] > gformpercent[i,9], "decrease2100",
                                      ifelse(gformpercent[i,6] == gformpercent[i,9], "stable", NA)))
  
}

#trend s3
for (i in 1:nrow(gformpercent)){
  gformpercent[i,16] <- ifelse(gformpercent[i,10] < gformpercent[i,13], "increase2100",
                               ifelse(gformpercent[i,10] > gformpercent[i,13], "decrease2100",
                                      ifelse(gformpercent[i,10] == gformpercent[i,13], "stable", NA)))
  
}

colnames(gformpercent) [14:16] <- c("trend_s1", "trend_s2", "trend_s3")

summary(as.factor(gformpercent$trend_s1))
summary(as.factor(gformpercent$trend_s2))
summary(as.factor(gformpercent$trend_s3))
