#p3ABC


#### packages ####
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(tidyverse)
library(ggplot2)

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{data}/predArray.rda")) #array
load(glue::glue("{data}/grid_25km.rda")) #grid
load(glue::glue("{data}/spTable.rda")) #spTable

spArray <- predArray[,,,,2]

gridsf <- st_as_sf(grid)

com1 <- gridsf %>%
  mutate(sp = apply(spArray[,,1,1], 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))


summary(com1$sp)
        
summary(apply(predArray[,,1,1,1],2, sum, na.rm = TRUE))

pres <- as.data.frame(predArray[,,1,1,1])

com2 <- gridsf %>%
  mutate(sp = apply(spArray[,,4,3], 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
summary(apply(spArray[,,4,3], 2, sum, na.rm = TRUE))

#unconstrained
spArray <- predArray[,,,,1]

com3 <- gridsf %>%
  mutate(sp = apply(spArray[,,4,3], 2, sum, na.rm = TRUE)) %>%
  dplyr::select(sp) %>%
  st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
summary(apply(spArray[,,4,3], 2, sum, na.rm = TRUE))

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

print(p11)
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
  nrow = 1,
  rel_widths = c(1,1,1,0.2)  # Adjust relative widths to minimize space further
)
# Print the final combined plot with a compact legend
print(combined_plot_with_legend)
ABC <- combined_plot_with_legend

#### boxplot ####

#data
summary(apply(predArray[,,1,1,1],2, sum, na.rm = TRUE))
sp1 <- apply(predArray[,,1,1,2],2, sum, na.rm = TRUE)
sp2 <- apply(predArray[,,4,3,2],2, sum, na.rm = TRUE)
sp3 <- apply(predArray[,,4,3,1],2, sum, na.rm = TRUE)

spdata <- as.data.frame(cbind(sp1,sp2,sp3))
spdata$cell <- 1:38657

boxplot(data=spdata[,1:3], x=spdata$cell)

str(spdata)

library(ggplot2)
library(tidyr)


# Reshape the data from wide to long format
spdata_long <- pivot_longer(spdata, cols = c(sp1, sp2, sp3), names_to = "sp", values_to = "value")


spdata_long <- spdata_long %>%
  mutate(sp = case_when(
    sp == "sp1" ~ "2010",
    sp == "sp2" ~ "2100, cons",
    sp == "sp3" ~ "2100, uncons",
    TRUE ~ as.character(sp)  # Keep the original value as a string if none of the above conditions match
  ))

# Create the boxplot
dataplot <- ggplot(spdata_long, aes(x = sp, y = value, fill = sp)) +
  geom_boxplot() +
  labs(x = "year and scenario", y = "species", title = "D) richness", fill = "data") +
  theme_minimal() +
  scale_fill_manual(values = c("2010" = "#66B2FF", "2100, uncons" = "#66CC66", "2100, cons" = "#FF9999"))+
  theme(      plot.title = element_text(size = 16, hjust = 0.5),
           #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
           axis.text.x = element_blank(),
           axis.title.x = element_text(size = 14),
           axis.text.y = element_text(size = 14),
           axis.title.y = element_text(size = 14))


print(dataplot)

summary(spdata$sp1)
summary(spdata$sp2)
summary(spdata$sp3)
