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
path <- "bioing/user/roschw"
path <- "C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/timeline_viridis"
#load(glue::glue("{path}/cormatrix_cells.rda"))
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

load(glue::glue("{data}/sigma1Array.rda")) #array
sp1 <- as.data.frame(spArray[,,1,1])
sp1 <- as.data.frame(t(sp1))

membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))



#### ndist ####

data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"
load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid
grid <- as.data.frame(spArrayMeta[[2]])

##### distance ####

center_point <- st_sfc(st_point(c(0, 0)), crs = st_crs(spArrayMeta[[2]]))

distg <- spArrayMeta[[2]] %>%
  mutate(dist_n = st_distance(geometry, center_point))

colnames(distg)[5] <- "distN"
#cr <- as.data.frame(matrix(nrow=6, ncol=38657))
distN <- as.data.frame(distg[,5])

membercldat1$dist <- distg$distN
membercldat4$dist <- distg$distN



# Load necessary libraries
library(ggplot2)
library(dplyr)
library(units)


# Assign categories based on the presence of keywords
membercldat4$Category[grepl("tundra", membercldat4$new_group)] <- "tundra"
membercldat4$Category[grepl("taiga, eurasian", membercldat4$new_group)] <- "taiga, eurasian"
membercldat4$Category[grepl("taiga, nearctic", membercldat4$new_group)] <- "taiga, nearctic"

membercldat1$Category = membercldat1$new_group


#membercldat4$Category[grepl("tundra", membercldat4$new_group)] <- "tundra"
#membercldat4$Category[grepl("taiga", membercldat4$new_group)] <- "taiga"
#membercldat4$Category[grepl("taiga", membercldat4$new_group)] <- "taiga, nearctic"

#membercldat1$Category = membercldat1$new_group
# membercldat1$Category[grepl("tundra", membercldat1$new_group)] <- "tundra"
# membercldat1$Category[grepl("taiga", membercldat1$new_group)] <- "taiga"


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

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# # Calculate summary statistics for year1
# year1_summary <- m1 %>%
#   group_by(Category) %>%
#   summarize(
#     min = min(dist),
#     Q1 = quantile(dist, 0.25),
#     median = median(dist),
#     Q3 = quantile(dist, 0.75),
#     max = max(dist),
#     .groups = 'drop'
#   ) %>%
#   mutate(year = "year1")
# 
# # Calculate summary statistics for year4  
# year4_summary <- m4 %>%
#   group_by(Category) %>%
#   summarize(
#     min = min(dist),
#     Q1 = quantile(dist, 0.25),
#     median = median(dist),
#     Q3 = quantile(dist, 0.75),
#     max = max(dist),
#     .groups = 'drop'
#   ) %>%
#   mutate(year = "year4")
# 
# # Combine the summaries into one data frame
# summary_combined <- bind_rows(year1_summary, year4_summary)
# 
# library(dplyr)
# 
# # Assuming summary_combined is your tibble
# # Step 1: Filter data for year1 and year3
# year1_data <- summary_combined %>% filter(year == "year1")
# year4_data <- summary_combined %>% filter(year == "year4")
# 
# # Step 2: Join the data on Category
# joined_data <- year1_data %>%
#   inner_join(year4_data, by = "Category", suffix = c("_year1", "_year4"))
# 
# # Step 3: Calculate differences
# result <- joined_data %>%
#   mutate(
#     min_diff = min_year1 - min_year4,
#     Q1_diff = Q1_year1 - Q1_year4,
#     median_diff = median_year1 - median_year4,
#     Q3_diff = Q3_year1 - Q3_year4,
#     max_diff = max_year1 - max_year4
#   ) %>%
#   select(Category, min_diff, Q1_diff, median_diff, Q3_diff, max_diff)
# 
# # Step 4: Sort by Category in ascending order (default behavior)
# final_result <- result %>%
#   arrange(Category)
# 
# dat <- as.data.frame(t(final_result))
# colnames(dat) <- dat[1,]
# dat <- dat[-1,]
# str(dat)
# 
# # Reshape the data from wide to long format
# long_df <- dat %>%
#   pivot_longer(cols = everything(), names_to = "Category", values_to = "Values")
# 
# long_df2 <- long_df %>%
#   filter(Category == "tundra")
# 
# # Step 1: Remove the units
# long_df$Values <- gsub("\\s*\\[m\\]", "", long_df$Values)  # Remove the '[m]' part
# 
# # Step 2: Convert to numeric
# long_df$Values  <- as.numeric(long_df$Values)

m4$ddist = m4$dist - median(m1$dist)
summary(m1$dist)

m4a <- m4 %>%
  filter(Category != "tundra")

p1 <- ggplot(data = m4, aes(x = Category, y = ddist*0.001)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.size = 1, outlier.shape = 16) +
  labs(x = "community", y = "distance [km]", title = "sig1 all cells") +
  theme_minimal() +
  ylim(0,150) +
  theme(axis.text.x = element_blank())  # Remove y-axis title
 # element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

print(p1)

#summary(long_df)

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

#dif2 <- as.data.frame(rbind(taiga_eu_diff, taiga_ne_diff, tundra_diff))
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

p1 <- ggplot(data = m4, aes(x = Category, y = ddist * 0.001, fill = Category)) +
  geom_boxplot(color = "black", outlier.size = 1, outlier.shape = 16) +
  labs(x = "", y = "distance [km]") +
  theme_minimal() +
  labs(title = "E) 2100 constrained")+ 
  #labs(title = "B) all cells") +  # Update labels
  theme(
      axis.text.x = element_blank(), #element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Increase x-axis text size
      axis.text.y = element_text(size = 14),  # Increase y-axis text size
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(),
      plot.title = element_text(size = 16),  # Increase title
    #axis.text.y = element_blank(),  # Remove y-axis text
    #axis.ticks.y = element_blank(),  # Remove y-axis ticks
    #axis.title.y = element_blank(),   # Remove y-axis title
    #axis.title.x = element_blank(),
    #panel.grid.major = element_blank(),  # Remove major grid lines
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    #panel.background = element_rect(fill = "white"),  # Set white background
    #panel.border = element_blank(), # Remove the panel border
    legend.position = "none",
    axis.title.y = element_text(size = 14))+  # Remove the legend
  ylim(0, 150) +
  scale_fill_manual(values = c("tundra" = tundra, 
                               "taiga, nearctic" = taiga_ne, 
                               "taiga, eurasian" = taiga_eu)) # Map categories to colors

print(p1)


p2 <- ggplot(data = dif, aes(x = Category, y = delta * 0.001, fill = Category)) +
  geom_boxplot(color = "black", outlier.size = 1, outlier.shape = 16) +
  labs(x = "", y = "", title = "new cells") +
  theme_minimal() +
  labs(title = expression(paste("F) ", Delta, " 2100 constrained"))) + 
 # labs(title = "C) new cells") +  # Update labels
  theme(
    axis.text.x = element_blank(), #element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    plot.title = element_text(size = 16),  # Increase title
    #axis.text.y = element_blank(),  # Remove y-axis text
    #axis.ticks.y = element_blank(),  # Remove y-axis ticks
    #axis.title.y = element_blank(),   # Remove y-axis title
    #axis.text.x = element_blank(), 
    #axis.ticks.x = element_blank(), 
    #axis.title.x = element_blank(),
    #panel.grid.major = element_blank(),  # Remove major grid lines
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    #panel.background = element_rect(fill = "white"),  # Set white background
    #panel.border = element_blank(), # Remove the panel border
    legend.position = "none")+  # Remove the legend
  ylim(0, 150) +
  scale_fill_manual(values = c("taiga, nearctic" = taiga_ne, 
                               "taiga, eurasian" = taiga_eu)) 

print(p2)


eu <- subset(dif, dif$Category == "taiga, eurasian")
ne <- subset(dif, dif$Category == "taiga, nearctic")

shapiro.test(eu$delta)
t.test(eu$delta, ne$delta)
wilcox.test(eu$delta, ne$delta)

library(cowplot)
combined_plot <- plot_grid(p1, p2 + theme(legend.position="none"), ncol=2, nrow=1) 
print(combined_plot)

P <- combined_plot
P <- add_sub(P, "communities", hjust = 0.3)
bc <- ggdraw(P)
print(bc)
