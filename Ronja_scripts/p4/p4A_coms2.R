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
library(ggplot2)
library(glue)



path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{path}/grid_25km.rda")) #grid

gridsf <- st_as_sf(grid)

#### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[1] #tundra
tundra2 <- hcl.colors(6, "Greens 3")[2] #tundra
tundra3 <- hcl.colors(6, "Greens 3")[3] #tundra
tundra4 <- hcl.colors(6, "Greens 3")[4] #tundra

taiga_ne <- hcl.colors(6, "Purples")[1] #nearctic
taiga_ne2 <- hcl.colors(6, "Purples")[2] #nearctic
taiga_ne3 <- hcl.colors(6, "Purples")[3] #nearctic
taiga_ne4 <- hcl.colors(6, "Purples")[4] #nearctic

taiga_eu <- hcl.colors(6, "Blues 3")[1] #europe
taiga_eu_W <- hcl.colors(6, "Blues 3")[2] #europe
taiga_eu_E <- hcl.colors(6, "Blues 3")[3] #europe
taiga_eu_E2 <- hcl.colors(6, "Blues 3")[4] #europe


col_all <- c(taiga_eu, taiga_eu_E, taiga_eu_E2, taiga_eu_W,
             taiga_ne, taiga_ne2, taiga_ne3, taiga_ne4,
             tundra, tundra2, tundra3)
  
#### y4 ####
  ##### load new data ####
  membercl <- read.csv(glue::glue("{path}/membercl_cells_y4s3_sigma1.csv"), sep =",")
  unique(membercl$group)
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~  "taiga, nearctic'''",
      group == 6 ~ "tundra''",
      group == 13 ~ "taiga, eurasian, W",
      group == 17 ~ "taiga, nearctic''",
      group == 19 ~ "taiga, eurasian, E'",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  #write.csv(membercl, "membercl4.csv", row.names = F)
  
  col4 = c(taiga_eu_E2, taiga_eu_W, taiga_ne3, taiga_ne4, tundra3)
  
  #### join ####
  library(dplyr)
  library(sf)
  library(stars)
 
  sp_with_groups <- gridsf %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com4 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
  
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

  ##### load new data ####
  membercl <- read.csv(glue::glue("{path}/membercl_cells_y3s3_sigma1.csv"), sep =",")
  unique(membercl$group)
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ NA,
      group == 6 ~ "tundra'",
      group == 10 ~ "taiga, nearctic'",
      group == 11 ~ "taiga, eurasian, W",
      group == 19 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  #write.csv(membercl, "membercl3.csv", row.names = F)
  col3 = c(taiga_eu_E, taiga_eu_W, taiga_ne2, tundra2)
  
  ##### join ####

  sp_with_groups <- gridsf %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com3 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
  
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

  ##### load new data ####
  membercl <- read.csv(glue::glue("{path}/membercl_cells_y2s3_sigma1.csv"), sep =",")
  unique(membercl$group)
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ NA,
      group == 5 ~ "tundra",
      group == 10 ~ "taiga, nearctic",
      group == 11 ~ "taiga, eurasian",
      group == 19 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  #write.csv(membercl, "membercl2.csv", row.names = F)
  col2 = c(taiga_eu, taiga_eu_E, taiga_ne, tundra)
  com2 <- gridsf %>% mutate(sp = membercl$new_group) %>% dplyr::select(sp) %>%
    st_rasterize(.,st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_)) 
  
  unique(membercl$new_group)
  
  ##### join ####

  sp_with_groups <- gridsf %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com2 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
  
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
  
  ##### load new data ####
  membercl <- read.csv(glue::glue("{path}/membercl_cells_y1s3_sigma1.csv"), sep =",")
  
  
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
      group == 10 ~ NA,
      group == 12 ~ "taiga, eurasian",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  #write.csv(membercl, "membercl1.csv", row.names = F)
  col1 = c(taiga_eu, taiga_ne, tundra)
  
  ##### join ####

  sp_with_groups <- gridsf %>%
    left_join(membercl %>% select(label, new_group), by = join_by(Id == label))
  
  # Rasterizing the data with new_group
  com1 <- sp_with_groups %>%
    mutate(sp = as.factor(new_group)) %>%
    dplyr::select(sp) %>%
    st_rasterize(., st_as_stars(st_bbox(gridsf), dx = 25050, dy = 25050, values = NA_real_))
  
  
  ##### plotting ####
  plot1 <- ggplot() +
    geom_stars(data = com1, aes(fill = as.factor(sp))) +  # Use fill for raster colors
    scale_fill_manual(values = col1, na.value = "transparent", na.translate = FALSE) +  # Ignore NA values) +  # Set custom colors
    theme_minimal() +
    ggtitle("A)", subtitle = "2010")+
    labs(fill = "Cell Communities") +  # Update labels #title = "2010", 
    theme(
      plot.title = element_text(size = 16),
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

  
  #### combined ####
  library(cowplot)
  A <- plot_grid(plot1, plot2, plot3, plot4 , ncol = 4)
  print(A)
  

  
