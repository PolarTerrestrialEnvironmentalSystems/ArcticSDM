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

#### clustering ####

data <- "/home/pd/roschw001/R/data"

load("//smb.isipd.dmawi.de/projects/bioing/user/roschw/sdm_occurance_metadata_25km_dispersal_99.rda") #grid


#load grid
data <- "//smb.isipd.dmawi.de/projects/bioing/user/roschw/"
load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid
setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/extinction_lag/")

#### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[1] #tundra only 1 type
# tundra2 <- hcl.colors(6, "Greens 3")[2] #tundra
# tundra3 <- hcl.colors(6, "Greens 3")[3] #tundra
# tundra4 <- hcl.colors(6, "Greens 3")[4] #tundra

taiga_ne <- hcl.colors(6, "Purples")[1] #nearctic
taiga_ne2 <- hcl.colors(6, "Purples")[2] #nearctic
# taiga_ne3 <- hcl.colors(6, "Purples")[3] #nearctic
# taiga_ne4 <- hcl.colors(6, "Purples")[4] #nearctic

taiga_eu <- hcl.colors(6, "Blues 3")[1] #europe
taiga_eu_W <- hcl.colors(6, "Blues 3")[2] #europe
taiga_eu_E <- hcl.colors(6, "Blues 3")[3] #europe
# taiga_eu_E2 <- hcl.colors(6, "Blues 3")[4] #europe


col_all <- c(taiga_eu, taiga_eu_E, taiga_eu_E2, taiga_eu_W,
             taiga_ne, taiga_ne2, taiga_ne3, taiga_ne4,
             tundra, tundra2, tundra3)
  
#### y4 ####
  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/membercl_cells_y4s3_ne.csv"), sep =",")
  unique(membercl$group)
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic'", 
      group == 10 ~ "taiga, eurasian, W",
      group == 14 ~ "taiga, eurasian, E", #18
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
  write.csv(membercl, "membercl4.csv", row.names = F)
  
  col4 = c(taiga_eu_E, taiga_eu_W, taiga_ne2, tundra)
  
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

  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/membercl_cells_y3s3_ne.csv"), sep =",")
  unique(membercl$group)
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic'", 
      group == 10 ~ "taiga, eurasian, W",
      group == 14 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
 write.csv(membercl, "membercl3.csv", row.names = F)
  col3 = c(taiga_eu_E, taiga_eu_W, taiga_ne2, tundra)
  
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

  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/membercl_cells_y2s3_ne.csv"), sep =",")
  unique(membercl$group)
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "tundra",
      group == 9 ~ "taiga, nearctic'", 
      group == 10 ~ "taiga, eurasian, W",
      group == 18 ~ "taiga, eurasian, E",
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
 write.csv(membercl, "membercl2.csv", row.names = F)
  col2 = c(taiga_eu_E, taiga_eu_W, taiga_ne2, tundra)
  com2 <- spArrayMeta[[2]] %>% mutate(sp = membercl$new_group) %>% dplyr::select(sp) %>%
    st_rasterize(.,st_as_stars(st_bbox(spArrayMeta[[2]]), dx = 25050, dy = 25050, values = NA_real_)) 
  
  unique(membercl$new_group)
  
  ##### join ####
  getwd()
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
  
  ##### load new data ####
  membercl <- read.csv(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/cellclusters_timeline_sigma1/membercl_cells_y1s3_sigma1.csv"), sep =",")
  
  
  # Step 1: Identify groups with only one label
  single_label_groups <- names(which(table(membercl$group) == 1))
  
  # Step 2: Update the group values to 0 for these labels
  membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0
  
  unique(membercl$group)
  
  membercl <- membercl %>%
    mutate(new_group = case_when(
      group == 0 ~ NA,
      group == 1 ~ "taiga, nearctic",
      group == 6 ~ "tundra",
      group == 10 ~ NA,
      group == 11 ~ "taiga, eurasian",
      group == 31 ~ NA,
      TRUE ~ as.character(group)  # Keep the original value as a string if none of the above conditions match
    ))
 write.csv(membercl, "membercl1.csv", row.names = F)
  col1 = c(taiga_eu, taiga_ne, tundra)
  
  ##### join ####
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
    ggtitle("constrained, no extinction", subtitle = "2010")+
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

  
  #### combined ####
  library(cowplot)
  combined_plot_res <- plot_grid(plot1, plot2, plot3, plot4 , ncol = 4)
  print(combined_plot_res)
  
  getwd()
 # setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/timeline_viridis")
  ggsave(combined_plot, file ="4plot.png", dpi = 900)
  

  
#### legend ####
  names <- as.data.frame(c(com1$sp, com2$sp, com3$sp, com4$sp))
  n <- unique(names)
  n <- n[2:12,]
  nso <- as.data.frame(sort(n$`c(com1$sp, com2$sp, com3$sp, com4$sp)`))
  str(nso)
  nso$col <- col_all
  
  cc <- com[1:11,]
  cc$col <- col_all
  cc$sp <- nso$`sort(n$/`c(com1$sp, com2$sp, com3$sp, com4$sp)/`)`
  
  
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
  