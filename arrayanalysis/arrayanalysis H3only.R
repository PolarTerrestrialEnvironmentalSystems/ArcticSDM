#1 with table 12 rows
#2 with table 4 rows, per scenario
#3 with species list and tab, 12 rows
#4 h3 with grid

#### input ####

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


#binarArray[species, cells, year, scenario]

#### initiate  ####
coords <- as.data.frame(st_coordinates(a$geometry, dims = c("x", "y"))) #extract coords

nam <- c("y1s1", "y2s1", "y3s1", "y4s1", 
           "y1s2","y2s2","y3s2", "y4s2",
           "y1s3","y2s3","y3s3","y4s3")

min_values <- as.data.frame(matrix(nrow=nrow(spTable), ncol = 12))  # Initialize a vector to store maximum values
max_values <- as.data.frame(matrix(nrow=nrow(spTable), ncol = 12))  # Initialize a vector to store maximum values
median_values <- as.data.frame(matrix(nrow=nrow(spTable), ncol = 12))  # Initialize a vector to store maximum values
mean_values <- as.data.frame(matrix(nrow=nrow(spTable), ncol = 12))  # Initialize a vector to store maximum values
colnames(min_values) <- nam
colnames(max_values) <- nam
colnames(median_values) <- nam
colnames(mean_values) <- nam

#### calculation loop ####
for (sp in 1:nrow(spTable)) {
  
species <- as.data.frame(binarArray[sp,,,]) #per species
a <- as.data.frame(cbind(grid, species))
xysp <- as.data.frame(cbind(coords, species))


for (i in 1:12) {
  v = i + 2
# Filter the data frame to keep only rows where a is equal to 1
filtered_data <- subset(xysp, xysp[, v] == 1)

# Find the maximum value of x in the filtered data
min_values[sp,i] <- min(filtered_data$Y)
max_values[sp,i] <- max(filtered_data$Y)
median_values[sp,i] <- median(filtered_data$Y)
mean_values[sp,i] <- mean(filtered_data$Y)

}
}

#### deltas ####
max_values$deltamax <- max_values$y4s3 - max_values$y1s1
max_values$changemax <- ifelse(max_values$deltamax > 0, "N",
                            ifelse(max_values$deltamax < 0, "S",
                                   ifelse(max_values$deltamax == 0, "stable", NA)))

min_values$deltamin <- min_values$y4s3 - min_values$y1s1
min_values$changemin <- ifelse(min_values$deltamin > 0, "N",
                            ifelse(min_values$deltamin < 0, "S",
                                   ifelse(min_values$deltamin == 0, "stable", NA)))

mean_values$deltamean <- mean_values$y4s3 - mean_values$y1s1
mean_values$changemean <- ifelse(mean_values$deltamean > 0, "N",
                            ifelse(mean_values$deltamean < 0, "S",
                                   ifelse(mean_values$deltamean == 0, "stable", NA)))

median_values$deltamedian <- median_values$y4s3 - median_values$y1s1
median_values$changemean <- ifelse(median_values$deltamedian > 0, "N",
                            ifelse(median_values$deltamedian < 0, "S",
                                   ifelse(median_values$deltamedian == 0, "stable", NA)))


#### result table ####
change <- as.data.frame(cbind(max_values$changemax, min_values$changemin, max_values$deltamax, min_values$deltamin))
colnames(change) <- c("change.max", "change.min", "delta.max", "delta.min")

change$shift <- ifelse(change$delta.max == 0 & change$delta.min == 0, "none", 
                       ifelse(change$change.max == "S" & change$change.min == "N", "shrink",
                              ifelse(change$change.max == "stable" & change$change.min == "S", "S_expansion",
                                     ifelse(change$change.max == "N" & change$change.min == "stable", "N_expansion",
                                            ifelse(change$change.max == "N" & change$change.min == "N", "N_shift", 
                                                   ifelse(change$change.max == "S" & change$change.min == "S", "S_expansion",
                                                          ifelse(change$change.max == "S" & change$change.min == "stable", "Sw_shrink",
                                                                 ifelse(change$change.max == "stable" & change$change.min == "N", "Nw_shrink",NA))))))))  
change$species = spTable$species



#### H4 change of communities ####

memberships <- read.csv("membercl_y1s1.csv", sep=",")
colnames(memberships) <- c("community", "species")

change$community <- memberships$community

change  %>%
  filter(community == "1")
