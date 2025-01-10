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
library(units)

path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{path}/grid_25km.rda")) #grid
gridsf <- st_as_sf(grid)

load(glue::glue("{path}/predArray.rda")) #array
spArray <- predArray[,,,,2]
sp1 <- as.data.frame(spArray[,,1,1])
sp1 <- as.data.frame(t(sp1))

membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))

#### ndist ####


##### distance ####

center_point <- st_sfc(st_point(c(0, 0)), crs = st_crs(grid))

distg <- grid %>%
  mutate(dist_n = st_distance(geometry, center_point))

colnames(distg)[5] <- "distN"

distN <- as.data.frame(distg[,5])

membercldat1$dist <- distg$distN
membercldat4$dist <- distg$distN

# Assign categories based on the presence of keywords
membercldat4$Category[grepl("tundra", membercldat4$new_group)] <- "tundra"
membercldat4$Category[grepl("taiga, eurasian", membercldat4$new_group)] <- "taiga, eurasian"
membercldat4$Category[grepl("taiga, nearctic", membercldat4$new_group)] <- "taiga, nearctic"

membercldat1$Category = membercldat1$new_group

membercldat4$dist <- as.numeric(set_units(membercldat4$dist, NULL))
membercldat1$dist <- as.numeric(set_units(membercldat1$dist, NULL))

m1 = as.data.frame(cbind(membercldat1$label, membercldat1$Category, membercldat1$dist))
m1 = na.omit(m1)
colnames(m1) = c("label", "Category", "dist")
m1$dist = as.numeric(m1$dist)

m4 = as.data.frame(cbind(membercldat4$label, membercldat4$Category, membercldat4$dist))
m4 = na.omit(m4)
colnames(m4) = c("label", "Category", "dist")
m4$dist = as.numeric(m4$dist)

##### ndist delta ####

m4$ddist = m4$dist - median(m1$dist)
summary(m1$dist)

m4a <- m4 %>%
  filter(Category != "tundra")


##### new cells ####
##### simple ####

taiga1 <- m1 %>%
  filter(Category %in% "taiga")

taiga2 <- m4 %>%
  filter(Category %in% "taiga")

taiga_diff <- anti_join(taiga2, taiga1, by = "label")

taiga_diff$delta <- -taiga_diff$dist + median(taiga1$dist)

taiga_eu1 <- m1 %>%
  filter(Category %in% "taiga, eurasian")


tundra1 <- m1 %>%
  filter(Category %in% "tundra")

tundra2 <- m4 %>%
  filter(Category %in% "tundra")


tundra_diff <- anti_join(tundra2, tundra1, by = "label")
#category4 has the category with new cells

tundra_diff$delta <- -tundra_diff$dist + median(tundra1$dist)


dif <- as.data.frame(rbind(taiga_diff, tundra_diff))
##############################################

####### 3 coms #####
taiga_ne1 <- m1 %>%
  filter(Category %in% "taiga, nearctic")

taiga_ne2 <- m4 %>%
  filter(Category %in% "taiga, nearctic")

taiga_ne_diff <- anti_join(taiga_ne2, taiga_ne1, by = "label")

taiga_ne_diff$delta <- -taiga_ne_diff$dist + median(taiga_ne1$dist)

taiga_eu1 <- m1 %>%
  filter(Category %in% "taiga, eurasian")

taiga_eu2 <- m4 %>%
  filter(Category %in% "taiga, eurasian")

taiga_eu_diff <- anti_join(taiga_eu2, taiga_eu1, by = "label")

taiga_eu_diff$delta <- -taiga_eu_diff$dist + median(taiga_eu1$dist)

tundra1 <- m1 %>%
  filter(Category %in% "tundra")

tundra2 <- m4 %>%
  filter(Category %in% "tundra")


tundra_diff <- anti_join(tundra2, tundra1, by = "label")
#category4 has the category with new cells

tundra_diff$delta <- -tundra_diff$dist + median(tundra1$dist)

dif <- as.data.frame(rbind(taiga_eu_diff, taiga_ne_diff, tundra_diff))
dif <- as.data.frame(rbind(taiga_eu_diff, taiga_ne_diff))



#### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[2] #tundra
taiga_ne <- hcl.colors(6, "Purples")[2] #america
taiga_eu <- hcl.colors(6, "Blues 3")[2] #europe

C <- ggplot(data = m4, aes(x = Category, y = ddist * 0.001, fill = Category)) +
  geom_boxplot(color = "black", outlier.size = 1, outlier.shape = 16) +
  labs(x = "", y = "distance [km]") +
  theme_minimal() +
  labs(title = "C) 2100 constrained")+ 
  #labs(title = "B) all cells") +  # Update labels
  theme(
      axis.text.x = element_blank(), #element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Increase x-axis text size
      axis.text.y = element_text(size = 14),  # Increase y-axis text size
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      plot.title = element_text(size = 16),  # Increase title
    axis.title.y = element_text(size = 14),  # Remove the legend
  legend.position = "none")+
  ylim(0, 150) +
  scale_fill_manual(values = c("tundra" = tundra, 
                               "taiga, nearctic" = taiga_ne, 
                               "taiga, eurasian" = taiga_eu)) # Map categories to colors

print(C)


B <- ggplot(data = dif, aes(x = Category, y = delta * 0.001, fill = Category)) +
  geom_boxplot(color = "black", outlier.size = 1, outlier.shape = 16) +
  labs(x = "",  y = "migration distance [km]", title = "new cells") +
  theme_minimal() +
  labs(title = "B)") + 
 # labs(title = ") new cells") +  # Update labels
  theme(
    axis.text.x = element_blank(), #element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16),  # Increase title
    legend.position = "none")+  # Remove the legend
  ylim(0, 150) +
  scale_fill_manual(values = c("taiga, nearctic" = taiga_ne, 
                               "taiga, eurasian" = taiga_eu)) 

print(B)


eu <- subset(dif, dif$Category == "taiga, eurasian")
ne <- subset(dif, dif$Category == "taiga, nearctic")


#Anderson-Darling Test:
library(nortest)
ad_test_result <- ad.test(ne$delta)
print(ad_test_result)

wilcox.test(eu$delta, ne$delta)

