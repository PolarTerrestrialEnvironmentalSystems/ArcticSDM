#arrayanalysis H2

#source louvain_analysis
#https://stackoverflow.com/questions/49834827/louvain-community-detection-in-r-using-igraph-format-of-edges-and-vertices

#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)

library(psych)
library(igraph)

#### clustering ####

#data <- "/bioing/R/data"
data <- "/home/pd/roschw001/R/data"
#data <- "/bioing/user/roschw/"

load("//smb.isipd.dmawi.de/projects/bioing/user/roschw/sdm_occurance_metadata_25km_dispersal_99.rda") #grid


#load grid
data <- "//smb.isipd.dmawi.de/projects/bioing/user/roschw/"
load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid
grid <- as.data.frame(spArrayMeta[[2]])

library(ggplot2)
library(glue)


#### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[1] #tundra
tundra2 <- hcl.colors(6, "Greens 3")[2] #tundra
tundra3 <- hcl.colors(6, "Greens 3")[3] #tundra
tundra4 <- hcl.colors(6, "Greens 3")[4] #tundra

taiga_ne <- hcl.colors(6, "Purples")[1] #america
taiga_in <- hcl.colors(6, "Purples")[2] #tundra
taiga_ne3 <- hcl.colors(6, "Purples")[3] #tundra
taiga_ne4 <- hcl.colors(6, "Purples")[4] #tundra

taiga_eu <- hcl.colors(6, "Blues 3")[1] #europe
taiga_out <- hcl.colors(6, "Blues 3")[3] #europe

taiga_eu_E <- hcl.colors(6, "Teal Grn")[2] #europe
taiga_eu_E2 <- hcl.colors(6, "Teal Grn")[3] #europe

show_col(viridis_pal()(20))

unclustered = "black"

tundra <- viridis_pal()(20)[16] #tundra
tundra2 <- viridis_pal()(20)[17] #tundra
tundra3 <- viridis_pal()(20)[18] #tundra
tundra4 <- viridis_pal()(20)[19] #tundra

taiga_ne <- viridis_pal()(20)[1] #america
taiga_ne2 <- viridis_pal()(20)[2] #america
taiga_ne3 <- viridis_pal()(20)[3] #america
taiga_ne4 <- viridis_pal()(20)[4] #america

taiga_eu <- viridis_pal()(20)[5] #europe
taiga_eu_N <- viridis_pal()(20)[6] #europe
taiga_eu_S <- viridis_pal()(20)[7] #europe
taiga_eu_S1 <- viridis_pal()(20)[8] #europe

other = "pink"


col_all <- c(taiga_eu, taiga_ne, tundra, taiga_eu_E, taiga_in,
             taiga_out, 
             tundra2, taiga_eu_E2)



col_sankey2 = c(other, taiga_eu, taiga_eu_N, taiga_eu_S, 
                taiga_ne, taiga_ne2, tundra, tundra2, 
                tundra3, unclustered, taiga_eu_S1,taiga_ne3, taiga_ne4, tundra4)

  y = 4
  
#### y4 ####
  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/cellclusters_sig5/membercl_cells_y4s3_sig5.csv"), sep =",")
  
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  unique(membercl$group)
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "taiga, inner", #ne
      group == 8 ~ "tundra'",
      group == 13 ~ "taiga, outer", #eu
      group == 14 ~ "taiga, eurasian, E'",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, "membercl4.csv", row.names = F)
  col4 = c(taiga_eu_E2, taiga_in, taiga_out, tundra2)
  
  
  #### join ####
  library(dplyr)
  library(sf)
  library(stars)
  
  # Assuming spArrayMeta[[2]] has a column named 'label' that matches membercl$label
  sp_with_groups <- spArrayMeta[[2]] %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com4 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))
  
 ##### plotting ####
  plot4 <- ggplot() +
    geom_stars(data = com4, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col4, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("", subtitle = "2100")+
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.title.y = element_blank(),   # Remove y-axis title
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.border = element_blank(), # Remove the panel border
      legend.position = "none"  # Remove the legend
    )
  print(plot4)
  
  #y3 ################################################################################
  y = 3
  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/cellclusters_sig5/membercl_cells_y3s3_sig5.csv"), sep =",")
  
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  unique(membercl$group)
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "taiga, inner", #ne
      group == 6 ~ "tundra'",
      group == 13 ~ "taiga, outer", #eu
      group == 14 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, "membercl3.csv", row.names = F)
  col3 = c(taiga_eu_E2, taiga_in, taiga_out, tundra2)
  
  ##### join ####
  # Assuming spArrayMeta[[2]] has a column named 'label' that matches membercl$label
  sp_with_groups <- spArrayMeta[[2]] %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com3 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))
  
  ##### plotting ####
  plot3 <- ggplot() +
    geom_stars(data = com3, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col3, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("", subtitle = "2070")+
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.title.y = element_blank(),   # Remove y-axis title
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.border = element_blank(), # Remove the panel border
      legend.position = "none"  # Remove the legend
    )
  print(plot3)
  
  
  #y2 ################################################################################
  y = 2
  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/cellclusters_sig5/membercl_cells_y2s3_sig5.csv"), sep =",")
  
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  unique(membercl$group)
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ NA,
      group == 6 ~ "tundra",
      group == 10 ~ "taiga, nearctic",
      group == 12 ~ "taiga, eurasian",
      group == 17 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, "membercl2.csv", row.names = F)
  col2 = c(taiga_eu, taiga_eu_E,taiga_ne, tundra)
  com2 <- spArrayMeta[[2]] %>% mutate(sp = membercl$new_group) %>% dplyr::select(sp) %>%
    st_rasterize(.,st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_)) 
  
  unique(membercl$new_group)
  ##### join ####
  
  # Assuming spArrayMeta[[2]] has a column named 'label' that matches membercl$label
  sp_with_groups <- spArrayMeta[[2]] %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com2 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))
  
  ##### plotting ####
  plot2 <- ggplot() +
    geom_stars(data = com2, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col2, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("", subtitle = "2040")+
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.title.y = element_blank(),   # Remove y-axis title
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.border = element_blank(), # Remove the panel border
      legend.position = "none"  # Remove the legend
    )
  print(plot2)
  #y1 ########################################################################
  
  y = 1
  #### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/cellclusters_sig5/membercl_cells_y1s3_sig5.csv"), sep =",")
  
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  unique(membercl$group)
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic",
      group == 10 ~ "taiga, eurasian",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, "membercl1.csv", row.names = F)
  col1 = c(taiga_eu, taiga_ne, tundra)
  
  #### join ####
  # Assuming spArrayMeta[[2]] has a column named 'label' that matches membercl$label
  sp_with_groups <- spArrayMeta[[2]] %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com1 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_))
  
  ##### plotting ####
  plot1 <- ggplot() +
    geom_stars(data = com1, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col1, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("B) timeline unconstrained", subtitle = "2010") +
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      axis.text.y = element_blank(),  # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      axis.title.y = element_blank(),   # Remove y-axis title
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.border = element_blank(), # Remove the panel border
      legend.position = "none"  # Remove the legend
    )
  print(plot1)

  
  ###########################
 
  
  #### combined ####
  library(cowplot)
  combined_plot_unres <- plot_grid(plot1, plot2, plot3, plot4 , ncol = 4)
  print(combined_plot_unres)
  
  getwd()
  setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/cellclusters_sig5")
  ggsave(combined_plot, file ="4plot.png", dpi = 900)
  
 
  
#### legend ####
  names <- as.data.frame(c(com1$sp, com2$sp, com3$sp, com4$sp))
  n <- unique(names)
  n <- n[2:12,]
  nso <- as.data.frame(sort(n$`c(com1$sp, com2$sp, com3$sp, com4$sp)`))
  
  nso <- nso[-5,]
  
  
  
  com <- grid %>%
    mutate(sp = as.factor(membercl$new_group)) %>% st_as_sf()
  
  cc <- com[1:8,]
  cc$col <- col_all
  cc$sp <- nso
  
  
  cc$sp <- c("taiga, eurasian", "taiga, eurasian, E", "taiga, eurasian, W", "taiga, eurasian, E'",   
            "taiga, nearctic", "taiga, nearctic'","taiga, nearctic''", "taiga, nearctic'''",         
            "tundra", "tundra'", "tundra''" )
  
  plotlegend <- ggplot(cc) +
    geom_sf(aes(color = as.factor(sp)), size = 0.001) +
    scale_color_manual(values = cc$col) +
    # Apply custom colors
    theme_minimal() +
    labs(title = glue::glue("y{y}s3 sig5 dispersal"), color = "cell communities")+
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  
  
  print(plotlegend)
  
  legend_plot <- ggplot(cc) +  # Use any dataset
    geom_sf(aes(color = as.factor(sp)), size = 0.001) +
    scale_color_manual(values = col_all) +  # Use the first color palette
    theme_minimal() +
    labs(color = "cell communities") +
    guides(color = guide_legend(override.aes = list(size = 5)))+
    theme(legend.position = "bottom")  # Position the legend
  
  
  # Create a separate plot for the legend using the 'cc' dataset
  legend_plot <- ggplot(cc) +  # Use your dataset for the legend
    geom_sf(aes(color = as.factor(sp)), size = 0.001) +
    scale_color_manual(values = col_all) +  # Use the color palette for the legend
    theme_minimal() +
    labs(color = "cell communities") +
    theme(legend.position = "bottom")  # Position the legend
  
  print(legend_plot)
  