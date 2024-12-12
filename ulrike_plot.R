
#plot ulrike

#### packages ####
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)



# Laden einer .tif-Datei als Vektor
bio1 <- rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/CHELSA data/CHELSA-Modern/CHELSA_Modern_Downsampling_5km/CHELSA-Modern/data/CHELSA_bio1_1981-2010_V.2.1.tif")
bio12 <- rast("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/environment/CHELSA data/CHELSA-Modern/CHELSA_Modern_Downsampling_5km/CHELSA-Modern/data/CHELSA_bio12_1981-2010_V.2.1.tif")

bio1$`CHELSA_bio1_1981-2010_V.2.1`
summary(bio1$`CHELSA_bio1_1981-2010_V.2.1`)

crs(bio1) = crs(spArrayMeta[[2]])
crs(bio1)

##########

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"
load(glue::glue("{data}/sigma1Array.rda")) #array
sig1 <- spArray
load(glue::glue("{data}/sig5Array.rda")) #array
sig5 <- spArray

load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid
grid <- as.data.frame(spArrayMeta[[2]])

##

sp = (apply(sig5[,,4,3], 2, sum) - apply(sig1[,,4,3], 2, sum))/apply(sig5[,,4,3], 2, sum)
sp1 = apply(sig1[,,4,3], 2, sum)
sp5 = apply(sig5[,,4,3], 2, sum)

# Extract the coordinates (points) from spArrayMeta[[2]]
coordinates <- st_coordinates(spArrayMeta[[2]])

# 3. Extract raster values at the spatial points
# Convert bio1$CHELSA_bio1_1981-2010_V.2.1 to a SpatRaster (if it isn't already)
raster_layer <- bio1$`CHELSA_bio1_1981-2010_V.2.1`
raster_layer12 <- bio12$`CHELSA_bio12_1981-2010_V.2.1`

# Extract values from the raster at the locations specified in the coordinates
extracted_values <- terra::extract(raster_layer, coordinates)
extracted_values12 <- terra::extract(raster_layer12, coordinates)

# 4. Add the extracted values to spArrayMeta[[2]]
spArrayMeta[[5]] <- spArrayMeta[[2]] %>%
  mutate(bio1_value = extracted_values)

spArrayMeta[[6]] <- spArrayMeta[[2]] %>%
  mutate(bio1_value = extracted_values12)

# View the updated spArrayMeta[[2]] with the extracted bio1 values
head(spArrayMeta[[5]])
head(spArrayMeta[[6]])

dat <- as.data.frame(matrix(nrow= 38657, ncol=1))
dat$V1 <- 1:38657
dat$bio1 <- unlist(spArrayMeta[[5]]$bio1_value)
dat$bio12 <-unlist(spArrayMeta[[6]]$bio1_value)
dat$ratio <- sp
dat$sig1 <- sp1
dat$sig5 <- sp5

colnames(dat)[1] <- "cell"

dat <- na.omit(dat)

#select subset of data
datq3 <- subset(dat, dat$sig1 < 20)

ggplot(data = datq3) + 
  geom_point(aes(x = bio1, y = bio12, color = sig1)) +
  scale_color_gradient(low = "blue", high = "red") +  # Adjust the color gradient as needed
  ylim(0,6000) +
  xlim(-25, 10)+
  labs(x = "annual mean temperature [°C]", y = "annual total precipitation [mm]", 
       title = "ratio < 0.2",
       color = "ratio") +   # Add labels for axes and color legend
  theme_minimal()  # Optional: Clean plot theme


ggplot(data = datq3) + 
  geom_point(aes(x = bio1, y = bio12, color = ratio)) +
  scale_color_gradient(low = "blue", high = "red",  limits = c(0, 1)) +  # Adjust the color gradient as needed
  ylim(0,6000) +
  xlim(-25, 10)+
  labs(x = "annual mean temperature [°C]", y = "annual total precipitation [mm]", 
       title = "ratio < 0.2",
       color = "ratio") +   # Add labels for axes and color legend
  theme_minimal()  # Optional: Clean plot theme

summary(dat$sig1)




ggplot(data = dat) + 
  geom_point(aes(x = bio1, y = bio12, color = sig5)) +
  scale_color_gradient(low = "blue", high = "red") +  # Adjust the color gradient as needed
  labs(x = "annual mean temperature [°C]", y = "annual total precipitation [mm]", 
       title = "unconstrained",
       color = "sig5") +   # Add labels for axes and color legend
  theme_minimal()  # Optional: Clean plot theme
