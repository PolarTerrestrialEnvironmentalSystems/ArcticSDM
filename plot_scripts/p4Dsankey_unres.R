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
library(mclust)

getwd()
path <- "C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/cellclusters_sig5"



membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))

membercldat1 = na.omit(membercldat1)
membercldat2 = na.omit(membercldat2)
membercldat3 = na.omit(membercldat3)
membercldat4 = na.omit(membercldat4)

# #### lookup
# 
# new_group <-  c(other, taiga_eu, c32, c33, c5b2, c5b3, taiga_ne, c5a1, c5a2, tundra, c22, c23, "black")
# new_group2 <- c("other", "taiga_eu", "taiga_eu_low", "taiga_eu_low'", "taiga_inner", "taiga_inner'", "taiga_ne",  
#               "taiga_outer", "taiga_outer'", "tundra", "tundra'", "tundra''", "unclustered")
# 
# 
# # Update membercldat1 with new_group2 based on new_group
# membercldat1 <- membercldat1 %>%
#   mutate(new_group2 = case_when(
#     new_group == "other" ~ "other",
#     new_group == "taiga_eu" ~ "taiga_eu",
#     new_group == "c32" ~ "taiga_eu_low", 
#     new_group == "c33" ~ "taiga_eu_low'",
#     new_group == "c5b2" ~ "taiga_inner", 
#     new_group == "c5b3" ~ "taiga_inner'",
#     new_group == "taiga_ne" ~ "taiga_ne",
#     new_group == "c5a1" ~ "taiga_outer",
#     new_group == "c5a2" ~ "taiga_outer'",
#     new_group == "tundra" ~ "tundra",
#     new_group == "c22" ~ "tundra'",
#     new_group == "c23" ~ "tundra''",
#     new_group == "unclustered" ~ "unclustered"
#   ))
# 
# setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/cellclusters_timeline_unres")
# write.csv(membercldat1, "membercl1n.csv", row.names = F)
# 
# 
# # Update membercldat1 with new_group2 based on new_group
# membercldat2 <- membercldat2 %>%
#   mutate(new_group2 = case_when(
#     new_group == "other" ~ "other",
#     new_group == "taiga_eu" ~ "taiga_eu",
#     new_group == "c32" ~ "taiga_eu_low", 
#     new_group == "c33" ~ "taiga_eu_low'",
#     new_group == "c5b2" ~ "taiga_inner", 
#     new_group == "c5b3" ~ "taiga_inner'",
#     new_group == "taiga_ne" ~ "taiga_ne",
#     new_group == "c5a1" ~ "taiga_outer",
#     new_group == "c5a2" ~ "taiga_outer'",
#     new_group == "tundra" ~ "tundra",
#     new_group == "c22" ~ "tundra'",
#     new_group == "c23" ~ "tundra''",
#     new_group == "unclustered" ~ "unclustered"
#   ))
# 
# setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/cellclusters_timeline_unres")
# write.csv(membercldat2, "membercl2n.csv", row.names = F)
# 
# 
# # Update membercldat1 with new_group2 based on new_group
# membercldat3 <- membercldat3 %>%
#   mutate(new_group2 = case_when(
#     new_group == "other" ~ "other",
#     new_group == "taiga_eu" ~ "taiga_eu",
#     new_group == "c32" ~ "taiga_eu_low", 
#     new_group == "c33" ~ "taiga_eu_low'",
#     new_group == "c5b2" ~ "taiga_inner", 
#     new_group == "c5b3" ~ "taiga_inner'",
#     new_group == "taiga_ne" ~ "taiga_ne",
#     new_group == "c5a1" ~ "taiga_outer",
#     new_group == "c5a2" ~ "taiga_outer'",
#     new_group == "tundra" ~ "tundra",
#     new_group == "c22" ~ "tundra'",
#     new_group == "c23" ~ "tundra''",
#     new_group == "unclustered" ~ "unclustered"
#   ))
# 
# setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/cellclusters_timeline_unres")
# write.csv(membercldat3, "membercl3n.csv", row.names = F)
# 
# 
# # Update membercldat1 with new_group2 based on new_group
# membercldat4 <- membercldat4 %>%
#   mutate(new_group2 = case_when(
#     new_group == "other" ~ "other",
#     new_group == "taiga_eu" ~ "taiga_eu",
#     new_group == "c32" ~ "taiga_eu_low", 
#     new_group == "c33" ~ "taiga_eu_low'",
#     new_group == "c5b2" ~ "taiga_inner", 
#     new_group == "c5b3" ~ "taiga_inner'",
#     new_group == "taiga_ne" ~ "taiga_ne",
#     new_group == "c5a1" ~ "taiga_outer",
#     new_group == "c5a2" ~ "taiga_outer'",
#     new_group == "tundra" ~ "tundra",
#     new_group == "c22" ~ "tundra'",
#     new_group == "c23" ~ "tundra''",
#     new_group == "unclustered" ~ "unclustered"
#   ))
# 
# setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/cellclusters_timeline_unres")
# write.csv(membercldat4, "membercl4n.csv", row.names = F)
# 
# 
# col_sankey <- c(c5b1, c5b2, c31, c32, c5a1, c5a2, c22, c21, unclustered, c5b3, c5a3, c23)
# 
# new_group2 <- c("other", "taiga_eu", "taiga_eu_low",  "taiga_inner",  "taiga_ne",  
#                 "taiga_outer", "taiga_outer'", "tundra", "tundra'", "taiga_eu_low'", "tundra''", "taiga_inner'", "unclustered")
# 
# col_sankey <-  c(other, taiga_eu, c32,  c5b2, taiga_ne, c5a1, c5a2, tundra, c22,  unclustered, c33,  c5b3, c23)
# 
# 
# year1_clusters <-  cbind(membercldat1$new_group2, membercldat1$label)
# year2_clusters <-  cbind(membercldat2$new_group2, membercldat2$label)
# year3_clusters <-  cbind(membercldat3$new_group2, membercldat3$label)
# year4_clusters <-  cbind(membercldat4$new_group2, membercldat4$label)

year1_clusters <-  cbind(membercldat1$new_group, membercldat1$label)
year2_clusters <-  cbind(membercldat2$new_group, membercldat2$label)
year3_clusters <-  cbind(membercldat3$new_group, membercldat3$label)
year4_clusters <-  cbind(membercldat4$new_group, membercldat4$label)

colnames(year1_clusters) <- c("cluster_group", "cell_number")
colnames(year2_clusters) <- c("cluster_group", "cell_number")
colnames(year3_clusters) <- c("cluster_group", "cell_number")
colnames(year4_clusters) <- c("cluster_group", "cell_number")

# Merge the data frames on cell_number
merged_data <- merge(year1_clusters, year2_clusters, by = "cell_number", suffixes = c("_year1", "_year2"))

merged_data <- merge(merged_data, year3_clusters, by = "cell_number")
colnames(merged_data)[4] <- "cluster_group_year3"
merged_data <- merge(merged_data, year4_clusters, by = "cell_number")
colnames(merged_data)[5] <- "cluster_group_year4"

colnames(merged_data) <- c("cell", "year1", "year2", "year3", "year4")

getwd()
write.csv(merged_data, "C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/clusters_dispersal.csv", row.names = F)

#### Sankey works ####
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
library(dplyr) # Also needed

merged_data <- read.csv("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/clusters_dispersal.csv")
colnames(merged_data) <- c("cells", "1", "2", "3", "4")

df2 <- merged_data %>%
  make_long("1", "2", "3", "4")


df2 <- df2 %>%
  mutate(year = case_when(
    x == 1 ~ 2010,
    x == 2 ~ 2040,
    x == 3 ~ 2070,
    x == 4 ~ 2100
    # Keep the original value as a string if none of the above conditions match
  ))

df2$year = as.factor(df2$year)

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

col_sankey2 <- c(taiga_eu, taiga_eu_E, taiga_in, taiga_ne,
             taiga_out, tundra,
             tundra2, taiga_eu_E2)

#########################################################
sankey_unres <- ggplot(df2, aes(x = as.factor(year), 
                          next_x = next_x, 
                          node = node, 
                          next_node = next_node,
                          fill = factor(node),
                          label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_manual(values = col_sankey2,
                    guide = guide_legend(title = "clusters"),
                    na.value = "transparent") +
  #theme_sankey(base_size = 16) +
  theme_sankey() +
  theme(legend.position = "none") +
  ggtitle("D) sankey unconstrained")+
  labs(x = "", y = "")


print(sankey_unres)
ggsave(sankey, file="sankey-plot.png", path ="//smb.isipd.dmawi.de/projects/bioing/user/roschw")

unique(df2$node)

#### combined ####

library(cowplot)
sankeys <- plot_grid(sankey_res, sankey_unres, ncol = 2)


combined <- plot_grid(combined_plot_res, combined_plot_unres, sankeys,nrow = 3)
                      

print(combined)
