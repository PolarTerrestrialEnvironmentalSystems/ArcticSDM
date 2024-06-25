#1 with table 12 rows
#2 with table 4 rows, per scenario
#3 with species list and tab, 12 rows
#4 h3 with grid

#### input ####


#meters

#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(tidyverse)
#input

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM"

load(glue::glue("{data}/sdm_arrays/spTable.rda"))
load(glue::glue("{data}/sdm_arrays/binaryArray_25km.rda"))
load(glue::glue("{data}/sdm_arrays/grid_25km.rda"))

load(glue::glue("{data}/sdm_occurance_array_array_25km.rda")) #only array

load(glue::glue("{data}/sdm_occurance_array_metadata_25km.rda")) #only meta
allsp <- as.data.frame(spArray[,,1,1])

a <- spArrayMeta[2]


distg <- grid %>%
  mutate(dist_n = st_distance(geometry, st_sfc(st_point(c(0, 0)), crs = st_crs(grid))))

plot(distg %>% dplyr::select(dist_n) %>% st_as_stars())


#grid
g <- st_as_stars(grid)
sppres <- as.data.frame(spArray[,,1,1])
sp = as.data.frame(spArrayList[[5]][,,1,1])

gs <- grid %>%
  mutate(values = sp)
gs  

richness <- grid %>% mutate(sp = apply(spArray[,,1,1], 2, sum))
                        
delta <- grid %>% mutate(a = apply(spArray[,,4,3], 2, sum),
                      b = apply(spArray[,,1,1], 2, sum),
                      c = b - a) %>% dplyr::select(c) %>% 
  st_rasterize(., st_as_stars(st_bbox(grid), dx = 25*1000, dy=25*1000, crs = st_crs(grid)))


plot(richness %>% dplyr::select(sp), cex=0.5, pch=16, main="species richness y1s1")

summary(richness)
plot(delta %>% dplyr::select(c), cex=0.5, pch=16)



#binarArray[species, cells, year, scenario]


#### initiate  ####
s <- as.data.frame(binarArray[1,,,]) 
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

#### calculation loop ####
for (sp in 1:nrow(spTable)) {
  
species <- as.data.frame(binarArray[sp,,,]) #per species
a <- as.data.frame(cbind(grid, species))
#Calculate decimal degree coordinates

X <- atan2(coords$Y, coords$X) * 180 / pi
Y <- atan2(coords$Y, coords$X * cos(X * pi/180)) * 180 / pi
xysp <- as.data.frame(cbind(X, Y, distg$distN, species))


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


change_avg <- as.data.frame(cbind(spTable$species, mean_values$changemean, median_values$changemean, mean_values$deltamean, median_values$deltamedian))
colnames(change_avg) <- c("species", "change.mean", "change.median", "delta.mean", "delta.median")

change_avg$summary <- ifelse(change_avg$change.mean == change_avg$change.median, change_avg$change.mean, "unclear")




#### H4 change of communities ####

memberships <- read.csv("membercl_y1s1.csv", sep=",")
colnames(memberships) <- c("community", "species")

changeTNS$community <- memberships$community

changeTNS  %>%
  filter(community == "3")
