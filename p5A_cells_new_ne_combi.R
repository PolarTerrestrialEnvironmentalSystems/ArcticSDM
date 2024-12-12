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
path <- "C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/cellclustering/extinction_lag"
     
#load(glue::glue("{path}/cormatrix_cells.rda"))
#data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

#load(glue::glue("{data}/sigma1Array.rda")) #array
#sp1 <- as.data.frame(spArray[,,1,1])
#sp1 <- as.data.frame(t(sp1))

membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))



summary(as.factor(membercldat4$new_group))


membercldat4 <- membercldat4 %>%
  mutate(new_group = case_when(
    new_group == "14" ~ "taiga, eurasian, E",
    TRUE ~ as.character(new_group)
  ))

membercldat1 <- membercldat1 %>%
  mutate(new_group = case_when(
    new_group == "other" ~ "NA",
    TRUE ~ as.character(new_group)
  ))


membercldat1 <-  na.omit(membercldat1)
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

# df1[1,1] <- "taiga, eurasian"
# 
# d1 <- as.data.frame(matrix(nrow=3, ncol =2))
# 
# d1$V1 <- c("tundra", "taiga, eurasian", "taiga, nearctic")
# d1$V2<- df1$Count1
# colnames(d1) = colnames(df1)

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

summary_table <- summary_table[-1,]
# Assuming percent_table is your data frame
summary_table_long_ne <- summary_table %>%
  pivot_longer(cols = starts_with("Count"), 
               names_to = "Count_Type", 
               values_to = "Value")

summary_table_long_ne$year <- c(2005, 2035, 2065, 2095, 2010, 2040, 2070, 2100, 2015, 2045, 2075, 2105)
summary_table_long_ne$year <- c(2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100,2010, 2040, 2070, 2100)


#### bar plot ####
a_ne <- ggplot(summary_table_long_ne, aes(x = Count_Type, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(taiga_eu, taiga_ne, tundra)) +
  labs(title = "G) Total cells, constrained (no extinction)",
       x = "",
       y = "colonized cells",
       fill = "Communities") +  # Change legend title to "Cluster"
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),  # Set legend title size
        plot.title = element_text(size = 16),    # Set plot title size
        axis.title.y = element_text(size = 14)) +
  scale_x_discrete(labels = c("Count1" = "2010", 
                              "Count2" = "2040", 
                              "Count3" = "2070", 
                              "Count4" = "2100"))

print(a_ne)

library(cowplot)
combined <- plot_grid(a, a_ne, ncol = 2, nrow=1)

eu <- subset(summary_table_long, summary_table_long$Category == "taiga, eurasian")
ne <- subset(summary_table_long, summary_table_long$Category == "taiga, nearctic")
tun <- subset(summary_table_long, summary_table_long$Category == "tundra")

shapiro.test(eu$Value)
wilcox.test(eu$ddist, tun$ddist)

#### combined ####
library(cowplot)
combined <- plot_grid(bc, a, ncol = 2, nrow=1)

print(combined)

