
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
library(geodata)

#########################
migp= "sigma1"

#### input ####
setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/hypotheses_results")
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid

path <- "//smb.isipd.dmawi.de/projects/bioing/user/slisovsk/ArcticSDM/Arrays_dispersalVector_Sigmoid"
files <- list.files(path, full.names = F, recursive = F)
spTable <- as.data.frame(files)
colnames(spTable) <- "species"

grid <- as.data.frame(spArrayMeta[[2]])

# Create a function to count the number of new 1 cells for each grid cell
count_new_ones <- function(present, future) {
  new_ones <- future - present
  new_ones[new_ones < 0] <- 0
  return(new_ones)
}


#### sig1 ####

load(glue::glue("{data}/sigma1Array.rda")) #array

present <- as.data.frame(spArray[, , 1, 1])
future <- as.data.frame(spArray[, , 4, 3])


new_array <- as.data.frame(count_new_ones(present, future))

n1 <- colSums(new_array)
plot(com1)

com1d <- spArrayMeta[[2]] %>%
  mutate(sp = n1) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))



##### plot ####

p1 <- ggplot() +
  geom_stars(data = com1d, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(option = "C", limits = c(0, 108), na.value = "transparent", name = "new species") +  # Apply Viridis color scale and handle NA values
  labs(title = expression(paste("D) ", Delta, " 2100 constrained"))) +  # Update labels    
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "right" , # Remove the legend
    options(repr.plot.width = 5, repr.plot.height = 4) 
  )
print(p1)
#max 108


#### sig5 ####

load(glue::glue("{data}/distconsArray.rda")) #array


present <- as.data.frame(spArray[, , 1, 1])
future <- as.data.frame(spArray[, , 4, 3])


new_array <- as.data.frame(count_new_ones(present, future))

n1 <- colSums(new_array)

summary(n1)

com5d <- spArrayMeta[[2]] %>%
  mutate(sp = n1) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))


##### plot ####

p2 <- ggplot() +
  geom_stars(data = com5d, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(option = "C", limits = c(0, 472), na.value = "transparent", name = "new species") +  # Apply Viridis color scale and handle NA values
  labs(title = expression(paste("E) ", Delta, " 2100 unconstrained"))) + # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    #plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "right", # Remove the legend +
    options(repr.plot.width = 5, repr.plot.height = 4) 
  )


print(p2)

#### both ####
library(cowplot)

# Combine plots with specified relative widths
DE <- plot_grid(
  NULL,   # Empty space for the first column (0.5)
  p1,     # Plot 1 in the second column (1)
  p2,     # Plot 2 in the third column (1)
  NULL,   # Empty space for the fourth column (0.5)
  ncol = 4,   # Specify number of columns
  rel_widths = c(0.4, 1.2, 1.2, 0.4) # Adjust relative widths
)

# Print the combined plot
print(DE)

ABCDE <- plot_grid(ABC, DE, nrow = 2)

print(ABCDE)

summary(com5d$sp)




