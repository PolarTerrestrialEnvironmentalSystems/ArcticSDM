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

path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))

membercldat1 = na.omit(membercldat1)
membercldat2 = na.omit(membercldat2)
membercldat3 = na.omit(membercldat3)
membercldat4 = na.omit(membercldat4)

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


#write.csv(merged_data, "clusters_dispersal.csv", row.names =F)

#### Sankey ####
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)


#### quick ####

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



col_sankey2 <- c(
             taiga_ne, taiga_ne2, taiga_eu, taiga_eu_E, taiga_eu_W,
             tundra, tundra2,taiga_ne4, taiga_ne3, taiga_eu_E2,  tundra3)
 
df2 <- df2 %>% mutate_all(~ gsub("taiga", "boreal forest", .))

df2 <- df2 %>% mutate_all(~ gsub("eurasian", "palearctic", .))



df2 <- df2 %>%
  mutate(label = ifelse(x == 1, node, NA))

unique(membercldat1$new_group)

C <- ggplot(df2, aes(x = as.factor(year), 
                          next_x = next_x, 
                          node = node, 
                          next_node = next_node,
                          fill = factor(node),
                          label = label)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_manual(values = col_sankey2,
                    guide = guide_legend(title = "clusters"),
                    na.value = "transparent") +
 theme_sankey() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16)) +
  ggtitle("C)")+
  labs(x = "", y = "") 



print(C)
