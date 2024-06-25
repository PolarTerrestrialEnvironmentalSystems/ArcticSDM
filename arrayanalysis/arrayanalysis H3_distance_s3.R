!#1 with table 12 rows
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
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

load(glue::glue("{data}/sdm_occurance_array_array_25km.rda")) #only array
load(glue::glue("{data}/sdm_occurance_array_metadata_25km.rda")) #only meta
load(glue::glue("{data}/sdm_arrays/grid_25km.rda")) #working grid

spTable <- as.data.frame(spArrayMeta[[1]])
colnames(spTable) <- "species"
#grid <- as.data.frame(spArrayMeta[[2]])



#spArray[species, cells, year, scenario]


#### initiate  ####
s <- as.data.frame(spArray[1,,,]) 
a <- as.data.frame(cbind(grid, s))
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

#### distance ####

distg <- grid %>%
  mutate(dist_n = st_distance(geometry, st_sfc(st_point(c(0, 0)), crs = st_crs(grid))))

plot(distg %>% dplyr::select(dist_n) %>% st_as_stars())


colnames(distg)[5] <- "distN"

#### calculation loop ####
for (sp in 1:nrow(spTable)) {
  
species <- as.data.frame(spArray[sp,,,]) #per species
a <- as.data.frame(cbind(grid, species))
#Calculate decimal degree coordinates

X <- atan2(coords$Y, coords$X) * 180 / pi
Y <- atan2(coords$Y, coords$X * cos(X * pi/180)) * 180 / pi
xysp <- as.data.frame(cbind(X, Y, distg$distN, species))
colnames(xysp)[3] <- "distN"

for (i in 1:12) {
  v = i + 3
# Filter the data frame to keep only rows where a is equal to 1
filtered_data <- subset(xysp, xysp[, v] == 1)

# Find the maximum value of x in the filtered data
min_values[sp,i] <- min(filtered_data$distN)
max_values[sp,i] <- max(filtered_data$distN)
median_values[sp,i] <- median(filtered_data$distN)
mean_values[sp,i] <- mean(filtered_data$distN)

}
}

#### deltas ####
##### max #####
max_values$deltamax1 <- max_values$y4s1 - max_values$y1s1
max_values$deltamax2 <- max_values$y4s2 - max_values$y1s1
max_values$deltamax3 <- max_values$y4s3 - max_values$y1s1
max_values$changemax1 <- ifelse(max_values$deltamax1 < 0, "N",
                            ifelse(max_values$deltamax1 > 0, "S",
                                   ifelse(max_values$deltamax1 == 0, "stable", NA)))
max_values$changemax2 <- ifelse(max_values$deltamax2 < 0, "N",
                               ifelse(max_values$deltamax2 > 0, "S",
                                      ifelse(max_values$deltamax2 == 0, "stable", NA)))
max_values$changemax3 <- ifelse(max_values$deltamax3 < 0, "N",
                               ifelse(max_values$deltamax3 > 0, "S",
                                      ifelse(max_values$deltamax3 == 0, "stable", NA)))

##### min #####
min_values$deltamin1 <- min_values$y4s1 - min_values$y1s1
min_values$deltamin2 <- min_values$y4s2 - min_values$y1s1
min_values$deltamin3 <- min_values$y4s3 - min_values$y1s1
min_values$changemin1 <- ifelse(min_values$deltamin1 < 0, "N",
                                ifelse(min_values$deltamin1 > 0, "S",
                                       ifelse(min_values$deltamin1 == 0, "stable", NA)))
min_values$changemin2 <- ifelse(min_values$deltamin2 < 0, "N",
                                ifelse(min_values$deltamin2 > 0, "S",
                                       ifelse(min_values$deltamin2 == 0, "stable", NA)))
min_values$changemin3 <- ifelse(min_values$deltamin3 < 0, "N",
                                ifelse(min_values$deltamin3 > 0, "S",
                                       ifelse(min_values$deltamin3 == 0, "stable", NA)))

##### mean #####
mean_values$deltamean1 <- mean_values$y4s1 - mean_values$y1s1
mean_values$deltamean2 <- mean_values$y4s2 - mean_values$y1s1
mean_values$deltamean3 <- mean_values$y4s3 - mean_values$y1s1
mean_values$changemean1 <- ifelse(mean_values$deltamean1 < 0, "N",
                                ifelse(mean_values$deltamean1 > 0, "S",
                                       ifelse(mean_values$deltamean1 == 0, "stable", NA)))
mean_values$changemean2 <- ifelse(mean_values$deltamean2 < 0, "N",
                                ifelse(mean_values$deltamean2 > 0, "S",
                                       ifelse(mean_values$deltamean2 == 0, "stable", NA)))
mean_values$changemean3 <- ifelse(mean_values$deltamean3 < 0, "N",
                                ifelse(mean_values$deltamean3 > 0, "S",
                                       ifelse(mean_values$deltamean3 == 0, "stable", NA)))

##### median #####
median_values$deltamedian1 <- median_values$y4s1 - median_values$y1s1
median_values$deltamedian2 <- median_values$y4s2 - median_values$y1s1
median_values$deltamedian3 <- median_values$y4s3 - median_values$y1s1
median_values$changemedian1 <- ifelse(median_values$deltamedian1 < 0, "N",
                                ifelse(median_values$deltamedian1 > 0, "S",
                                       ifelse(median_values$deltamedian1 == 0, "stable", NA)))
median_values$changemedian2 <- ifelse(median_values$deltamedian2 < 0, "N",
                                ifelse(median_values$deltamedian2 > 0, "S",
                                       ifelse(median_values$deltamedian2 == 0, "stable", NA)))
median_values$changemedian3 <- ifelse(median_values$deltamedian3 < 0, "N",
                                ifelse(median_values$deltamedian3 > 0, "S",
                                       ifelse(median_values$deltamedian3 == 0, "stable", NA)))

#### result table ####
change <- as.data.frame(cbind(max_values$changemax1,max_values$changemax2, max_values$changemax3, 
                              min_values$changemin1, min_values$changemin2, min_values$changemin3,
                              max_values$deltamax1, max_values$deltamax2, max_values$deltamax3,
                              min_values$deltamin1, min_values$deltamin2, min_values$deltamin3))

colnames(change) <- c("change.max1", "change.max2", "change.max3",
                      "change.min1", "change.min2", "change.min3",
                      "delta.max1", "delta.max2", "delta.max3",
                      "delta.min1", "delta.min2", "delta.min3")

change$shift1 <- ifelse(change$change.max1 == "stable" & change$change.min1 == "stable", "none", 
                       ifelse(change$change.max1 == "S" & change$change.min1 == "N", "shrink",
                              ifelse(change$change.max1 == "stable" & change$change.min1 == "S", "S_expansion",
                                     ifelse(change$change.max1 == "N" & change$change.min1 == "stable", "N_expansion",
                                            ifelse(change$change.max1 == "N" & change$change.min1 == "N", "N_shift", 
                                                   ifelse(change$change.max1 == "S" & change$change.min1 == "S", "S_shift",
                                                   ifelse(change$change.max1 == "N" & change$change.min1 == "S", "expansion",
                                                          ifelse(change$change.max1 == "S" & change$change.min1 == "stable", "Sw_shrink",
                                                                 ifelse(change$change.max1 == "stable" & change$change.min1 == "N", "Nw_shrink",NA))))))))) 

change$shift2 <- ifelse(change$delta.max2 == 0 & change$delta.min2 == 0, "none", 
                       ifelse(change$change.max2 == "S" & change$change.min2 == "N", "shrink",
                              ifelse(change$change.max2 == "stable" & change$change.min2 == "S", "S_expansion",
                                     ifelse(change$change.max2 == "N" & change$change.min2 == "stable", "N_expansion",
                                            ifelse(change$change.max2 == "N" & change$change.min2 == "N", "N_shift", 
                                                   ifelse(change$change.max2 == "S" & change$change.min2 == "S", "S_shift",
                                                   ifelse(change$change.max2 == "N" & change$change.min2 == "S", "expansion",
                                                          ifelse(change$change.max2 == "S" & change$change.min2 == "stable", "Sw_shrink",
                                                                 ifelse(change$change.max2 == "stable" & change$change.min2 == "N", "Nw_shrink",NA)))))))))  

change$shift3 <- ifelse(change$delta.max3 == 0 & change$delta.min3 == 0, "none", 
                       ifelse(change$change.max3 == "S" & change$change.min3 == "N", "shrink",
                              ifelse(change$change.max3 == "stable" & change$change.min3 == "S", "S_expansion",
                                     ifelse(change$change.max3 == "N" & change$change.min3 == "stable", "N_expansion",
                                            ifelse(change$change.max3 == "N" & change$change.min3 == "N", "N_shift", 
                                                   ifelse(change$change.max3 == "S" & change$change.min3 == "S", "S_shift", 
                                                   ifelse(change$change.max3 == "N" & change$change.min3 == "S", "expansion",
                                                          ifelse(change$change.max3 == "S" & change$change.min3 == "stable", "Sw_shrink",
                                                                 ifelse(change$change.max3 == "stable" & change$change.min3 == "N", "Nw_shrink",NA)))))))))  
change$species = spTable$species

summary(as.factor(change$shift1))
summary(as.factor(change$shift2))
summary(as.factor(change$shift3))

sum(174, 213,61, 69)

change_avg <- as.data.frame(cbind(spTable$species, 
                                  mean_values$changemean1, mean_values$changemean2, mean_values$changemean3,
                                  median_values$changemedian1, median_values$changemedian2, median_values$changemedian3,
                                  mean_values$deltamean1, mean_values$deltamean2, mean_values$deltamean3, 
                                  median_values$deltamedian1, median_values$deltamedian2, median_values$deltamedian3))

colnames(change_avg) <- c("species",
                          "change.mean1", "change.mean2", "change.mean3",
                          "change.median1", "change.median2", "change.median3",
                          "delta.mean",  "delta.mean",  "delta.mean",
                          "delta.median",  "delta.median",  "delta.median")

change_avg$summary1 <- ifelse(change_avg$change.mean1 == change_avg$change.median1, change_avg$change.mean1, "unclear")
change_avg$summary2 <- ifelse(change_avg$change.mean2 == change_avg$change.median2, change_avg$change.mean2, "unclear")
change_avg$summary3 <- ifelse(change_avg$change.mean3 == change_avg$change.median3, change_avg$change.mean3, "unclear")

summary(as.factor(change_avg$summary1))
summary(as.factor(change_avg$summary2))
summary(as.factor(change_avg$summary3))

write.csv2(change_avg, "H3_change_avg.csv", row.names = F)

plot(mean_values$y1s1/1000)

