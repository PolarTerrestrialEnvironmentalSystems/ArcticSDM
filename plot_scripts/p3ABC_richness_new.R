#H2 richness maps


#### packages ####
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(tidyverse)
library(ggplot2)

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"
load(glue::glue("{data}/sdm_occurance_array_array_25km.rda")) #only array

migp <- "sigma1"
setwd( "C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/hypotheses_results")
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

load(glue::glue("{data}/sigma1Array.rda")) #array

load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid

path <- "//smb.isipd.dmawi.de/projects/bioing/user/slisovsk/ArcticSDM/Arrays_dispersalVector_Sigmoid"
files <- list.files(path, full.names = F, recursive = F)
spTable <- as.data.frame(files)

colnames(spTable) <- "species"
grid <- as.data.frame(spArrayMeta[[2]])

com1 <- spArrayMeta[[2]] %>%
  mutate(sp = apply(spArray[,,1,1], 2, sum)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))

summary(apply(spArray[,,1,1], 2, sum))


com2 <- spArrayMeta[[2]] %>%
  mutate(sp = apply(spArray[,,4,3], 2, sum)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))
summary(apply(spArray[,,4,3], 2, sum))

load(glue::glue("{data}/sig5Array.rda")) #array

com3 <- spArrayMeta[[2]] %>%
  mutate(sp = apply(spArray[,,4,3], 2, sum)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))
summary(apply(spArray[,,4,3], 2, sum))

summary(com1$sp)
plot(com1)

sp = apply(spArray[,,4,3], 2, sum)

summary(sp)

##### plotting ####

  
library(ggplot2)
library(viridis)
library(stars)

p11 <- ggplot() +
  geom_stars(data = com1, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(limits = c(0, 853), na.value = "transparent", name = "species") +  # Apply Viridis color scale and handle NA values
  labs(title = "A) 2010") +  # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )

p2s1 <- ggplot() +
  geom_stars(data = com2, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(limits = c(0, 853), na.value = "transparent", name = "species") +
  labs(title = "B) 2100 constrained") +  # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )

p2s5 <- ggplot() +
  geom_stars(data = com3, aes(fill = sp)) +  # Use fill for raster colors
  scale_fill_viridis_c(limits = c(0, 853), na.value = "transparent", name = "species") +
  labs(title = "C) 2100 unconstrained") +  # Update labels
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white"),  # Set white background
    panel.border = element_blank(), # Remove the panel border
    legend.position = "none"  # Remove the legend
  )


summary(com2)

s1 = as.data.frame(spArray[,,1,1])
s1 = apply(spArray[,,1,1], 2, sum)
s2 = apply(spArray[,,4,3], 2, sum)
summary(s2)



###########
library(cowplot)

# Create a combined plot without a separate legend
combined_plot <- plot_grid(
  p11 + theme(legend.position = "none"),  # Remove legend from p1
  p2s1 + theme(legend.position = "none"),  # Remove legend from p2s1
  p2s5 + theme(legend.position = "none"),  # Remove legend from p2s5
  ncol = 3, 
  nrow = 1
)

print(combined_plot)

max(com1$sp)
com1

# Adjust the legend to reduce space around it
legend <- get_legend(
  p11 + 
    theme(
      legend.position = "right", 
      legend.box.margin = margin(0, 0, 0, 0),  # Remove margin around the legend box
      legend.spacing.y = unit(0, "cm"),        # Reduce vertical spacing between legend items
      legend.spacing.x = unit(0, "cm")         # Reduce horizontal spacing between legend items
    )
)

# Add the adjusted legend directly to the right of the combined plot
combined_plot_with_legend <- plot_grid(
  combined_plot, 
  legend, 
  rel_widths = c(3, 0.2)  # Adjust relative widths to minimize space further
)

combined_plot_with_legend <- plot_grid(
  p11 + theme(legend.position = "none"),  # Remove legend from p1
  p2s1 + theme(legend.position = "none"),  # Remove legend from p2s1
  p2s5 + theme(legend.position = "none"),  # Remove legend from p2s5
  legend, 
  ncol = 4,
  rel_widths = c(1,1,1,0.2)  # Adjust relative widths to minimize space further
)

combined_plot_with_legend <- plot_grid(
  p11 + theme(legend.position = "none"),  # Remove legend from p1
  p2s1 + theme(legend.position = "none"),  # Remove legend from p2s1
  p2s5,
  #p2s5 + theme(legend.position = "none"),  # Remove legend from p2s5
  legend, 
 NULL,
 NULL,
  p1,
    p2,
  ncol = 4,
  nrow=2)
  #rel_widths = c(1,1,1,0.2)  # Adjust relative widths to minimize space further

# Print the final combined plot with a compact legend
print(combined_plot_with_legend)
ABC <- combined_plot_with_legend

