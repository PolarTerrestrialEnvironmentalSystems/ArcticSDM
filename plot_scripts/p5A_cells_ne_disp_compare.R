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
#data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

#load(glue::glue("{data}/sigma1Array.rda")) #array
#sp1 <- as.data.frame(spArray[,,1,1])
#sp1 <- as.data.frame(t(sp1))

membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))


#### summary table ####

# Convert tables to data frames
df1 <- as.data.frame(table(membercldat1$new_group))
df2 <- as.data.frame(table(membercldat2$new_group))
df3 <- as.data.frame(table(membercldat3$new_group))
df4 <- as.data.frame(table(membercldat4$new_group))

# Set the names of the data frames to the first column (group names)
names(df1) <- c("Group", "Count1")
names(df2) <- c("Group", "Count2")
names(df3) <- c("Group", "Count3")
names(df4) <- c("Group", "Count4")

# Merge all data frames by Group
merged_df <- Reduce(function(x, y) merge(x, y, by = "Group", all = TRUE), list(df1, df2, df3, df4))

# Replace NA with NA (if needed) and display the final merged table
merged_df[is.na(merged_df)] <- NA

# Print the final table
print(merged_df)

# Sort the merged_df by the first column (Group) alphabetically
sorted_df <- merged_df[order(merged_df$Group, na.last = TRUE), ]

sorted_df$Category <- NA  # Initialize the new column with NA

# Assign categories based on the presence of keywords
sorted_df$Category[grepl("tundra", sorted_df$Group)] <- "tundra"
sorted_df$Category[grepl("taiga, eurasian", sorted_df$Group)] <- "Palearctic boreal forest"
sorted_df$Category[grepl("taiga, nearctic", sorted_df$Group)] <- "Nearctic boreal forest"
sorted_df$Category[grepl("other", sorted_df$Group)] <- "other"

sorted_df[is.na(sorted_df)] <- 0
summary_table <- aggregate(cbind(Count1, Count2, Count3, Count4) ~ Category,
                           data = sorted_df,
                           FUN = sum)



#### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[2] #tundra
taiga_ne <- hcl.colors(6, "Purples")[2] #america
taiga_eu <- hcl.colors(6, "Blues 3")[2] #europe


# Assuming percent_table is your data frame
summary_table_long <- summary_table %>%
  pivot_longer(cols = starts_with("Count"), 
               names_to = "Count_Type", 
               values_to = "Value")

#summary_table_long$year <- c(2005, 2035, 2065, 2095, 2010, 2040, 2070, 2100, 2015, 2045, 2075, 2105)

summary_table_long$year <- c(2010, 2040, 2070, 2100,2010, 2040, 2070, 2100,2010, 2040, 2070, 2100)

custom_colors <- c(taiga_ne, taiga_eu, tundra)



#### bar plot ####
a <- ggplot() +
  geom_point(data = summary_table_long, aes(x = year, y = Value, col = Category), stat = "identity", position = position_dodge()) +
  geom_line(data = summary_table_long, aes(x = year, y = Value, col = Category), stat = "identity", linetype = "dotted", position = position_dodge()) +
  geom_point(data = summary_table_long_ne, aes(x = year, y = Value, col = Category), shape = 17, stat = "identity", position = position_dodge()) +
  geom_line(data = summary_table_long_ne, aes(x = year, y = Value, col = Category), stat = "identity", linetype = "dashed", position = position_dodge()) +
  scale_fill_manual(values = c(taiga_eu, taiga_ne, tundra)) +
  labs(title = "Total cells, constrained vs no extinction",
       x = "",
       y = "colonized cells",
       color = "Communities") +  # Change legend title to "Cluster"
  scale_color_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks=c(2010, 2040, 2070, 2100), 
                     labels=c("2010", "2040", "2070", "2100")) + # Custom x-axis breaks and labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),  # Set legend title size
        legend.direction = "vertical",
        plot.title = element_text(size = 16),    # Set plot title size
        axis.title.y = element_text(size = 14)) +
  scale_x_discrete(labels = c("Count1" = "2010", 
                              "Count2" = "2040", 
                              "Count3" = "2070", 
                              "Count4" = "2100"))


print(a)

summary(as.factor(membercldat4$new_group))





eu <- subset(summary_table_long, summary_table_long$Category == "taiga, eurasian")
ne <- subset(summary_table_long, summary_table_long$Category == "taiga, nearctic")
tun <- subset(summary_table_long, summary_table_long$Category == "tundra")

shapiro.test(eu$Value)
wilcox.test(eu$ddist, tun$ddist)

#### combined ####
library(cowplot)
combined <- plot_grid(bc, a, ncol = 2, nrow=1)

print(combined)


###

summary_table_long$Value_ne <- summary_table_long_ne$Value
summary_table_long$delta <- summary_table_long$Value - 10000
summary_table_long$delta_ne <- summary_table_long$Value_ne - 10000
summary_table_long$delta_ratio <- summary_table_long$delta_ne/summary_table_long$delta*100
summary_table_long$delta_ratio2 <- summary_table_long$delta/summary_table_long$delta_ne*100

ggplot() +
  geom_point(data = summary_table_long, aes(x = year, y = delta_ratio, col = Category), stat = "identity", position = position_dodge()) +
  geom_line(data = summary_table_long, aes(x = year, y = delta_ratio, col = Category), stat = "identity", linetype = "dotted", position = position_dodge())

head(summary_table_long)

mean_delta_ratio2 <- summary_table_long %>%
  group_by(Category) %>%
  summarise(mean_delta_ratio = mean(delta_ratio2, na.rm = TRUE))

print(mean_delta_ratio2)

d4 <- subset(summary_table_long, summary_table_long$Count_Type == "Count4" | summary_table_long$Count_Type == "Count1")

(d4[2,5] -d4[1,3])/ (d4[2,3] -d4[1,3])  #nearc value 4 - #nearc value 1 /nearc value_ne 4 - #nearc value 1
# 0.21 of increase nearc 

(d4[4,5] -d4[3,3])/ (d4[4,3] -d4[3,3])  #palearc value 4 - #palearc value 1 /palearc value_ne 4 - #palearc value 1
# 0.68 of increase palearc

(d4[6,5] -d4[5,3])/ (d4[6,3] -d4[5,3])  #tundra value 4 - #tundra value 1 /tundra value_ne 4 - #tundra value 1
# 0.5 of decrease tundra

