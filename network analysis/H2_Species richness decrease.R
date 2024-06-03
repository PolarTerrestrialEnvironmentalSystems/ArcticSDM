#Hypothesis 2: Species richness decrease

#### packages ####
library(stars)
library(tidyr)
library(terra)
library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

#### input ####

data <- "D:/Results"
#data <- "C:/Users/roschw001/Documents/R/SDM/SDM/festplatte_sub"

#speciesnames

flsl      <- list.files(data, 
                        recursive = FALSE, full.names = T)

speciesnames = as.data.frame(matrix(nrow = 1, ncol=1))

# Split the path by the "/" character
for (i in 1:length(flsl)){
  
  # Extract the last part of the path, which contains the species names
  path_parts <- unlist(strsplit(flsl[i], "/|_"))
  
  # Extract the last part of the path, which contains the species name
  speciesnames[i,1] <- path_parts[3]#10
  
}

colnames(speciesnames)[1] <- "species"

#### Present ####

##### binary #####
binaries <- list()

for(sp in 1:3){
#sp=2
species = speciesnames[sp,1]
  
#as stars (NA outside area also becomes 0 to enable futher calulcations)

binaries[[sp]] <- read_stars(glue::glue("{data}/{species}_MaxEnt_calibration.tif")) %>%
  setNames("present") %>%
  mutate(present = ifelse(present>0.5, 1, 0)) %>%
  mutate(present = ifelse(is.na(present), 0, present))

}

enbies <- list()
for(sp in 1:3){
species = speciesnames[sp,1]
enbies[[sp]] <- read_stars(glue::glue("{data}/{species}_MaxEnt_calibration.tif")) %>%
  setNames("present") %>%
  mutate(present = present)
}

##### richness plot #####
#stack species layers and sum (species richness plot)

#stack
stack <- do.call(c, binaries)
stack2 <- do.call(c, enbies)

st_apply(stack, c("present", "present.1"), sum)
stack
#sum
stack <- stack %>%
  mutate(sum = present + present.1 + present.2 + present.3 + present.4 + present.5 + present.6
         + present.7 + present.8 + present.9 + present.10 + present.11 + present.12 + present.13
         + present.14 + present.15 + present.16 + present.17 + present.18 + present.19)

#plot sum
plot(stack[21])
stack[21]

#### Future ####

load(glue::glue("{data}/Results/{species}/Predictions/pred_stars.rda")) #future

#for(year in 1:2) {
#  for(scen in 1:2) {
year = 3
scen = 2

fut1 <- pred_stars[scen][,,,year]

##### binary #####
binariesf <- list()

for(sp in 1:20){
  
  species = speciesnames[sp,1]
  
  binariesf[[sp]] <- read_stars(glue::glue("{data}/{species}_MaxEnt_calibration.tif")) %>%
    setNames("present") %>%
    mutate(future = ifelse(present>0.5, 1, 0)) %>%
    mutate(future = ifelse(is.na(present), 0, present))
}

##### richness plot #####

#stack
stackf <- do.call(c, binaries)
stackf2 <- do.call(c, enbies)

#sum
stackf <- stackf %>%
  mutate(sum = future + future.1 + future.2 + future.3 + future.4 + future.5 + future.6
         + future.7 + future.8 + future.9 + future.10 + future.11 + future.12 + future.13
         + future.14 + future.15 + future.16 + future.17 + future.18 + future.19)

#plot sum
plot(stackf[21])
stackf[21]


