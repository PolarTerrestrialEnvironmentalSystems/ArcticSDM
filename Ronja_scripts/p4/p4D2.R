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
     
#### dispersal normal ####
membercldat1 <- read.csv(glue::glue("{path}/membercl1.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4.csv"))


##### summary table ####

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


##### colours ####
library(scales)

tundra <- hcl.colors(6, "Greens 3")[2] #tundra
taiga_ne <- hcl.colors(6, "Purples")[2] #america
taiga_eu <- hcl.colors(6, "Blues 3")[2] #europe


# Assuming percent_table is your data frame
summary_table_long <- summary_table %>%
  pivot_longer(cols = starts_with("Count"), 
               names_to = "Count_Type", 
               values_to = "Value")

summary_table_long$year <- c(2010, 2040, 2070, 2100,2010, 2040, 2070, 2100,2010, 2040, 2070, 2100)

custom_colors <- c(taiga_ne, taiga_eu, tundra)

#### no extinction ####

membercldat1 <- read.csv(glue::glue("{path}/membercl1_ne.csv"))
membercldat2 <- read.csv(glue::glue("{path}/membercl2_ne.csv"))
membercldat3 <- read.csv(glue::glue("{path}/membercl3_ne.csv"))
membercldat4 <- read.csv(glue::glue("{path}/membercl4_ne.csv"))

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

summary(as.factor(membercldat1$new_group))
##### summary table ####

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

summary_table_long_ne <- summary_table %>%
  pivot_longer(cols = starts_with("Count"), 
               names_to = "Count_Type", 
               values_to = "Value")

summary_table_long_ne$year <- c(2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100)



#### comparison plot ####
D <- ggplot() +
  geom_point(data = summary_table_long, aes(x = year, y = Value, col = Category), stat = "identity", position = position_dodge()) +
  geom_line(data = summary_table_long, aes(x = year, y = Value, col = Category), stat = "identity", linetype = "dotted", position = position_dodge()) +
  geom_point(data = summary_table_long_ne, aes(x = year, y = Value, col = Category), shape = 17, stat = "identity", position = position_dodge()) +
  geom_line(data = summary_table_long_ne, aes(x = year, y = Value, col = Category), stat = "identity", linetype = "dashed", position = position_dodge()) +
  scale_fill_manual(values = c(taiga_eu, taiga_ne, tundra)) +
  labs(title = "D)",
       x = "",
       y = "colonized cells")+
      # color = "Communities") +  # Change legeend title to "Cluster"
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks=c(2010, 2040, 2070, 2100), 
                     labels=c("2010", "2040", "2070", "2100")) + # Custom x-axis breaks and labels
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16))

print(D)

#### estimations ####

summary_table_long$Value_ne <- summary_table_long_ne$Value

summary(as.factor(membercldat4$new_group))

d4 <- subset(summary_table_long, summary_table_long$Count_Type == "Count4" | summary_table_long$Count_Type == "Count1")

(d4[2,5] -d4[1,3])/ (d4[2,3] -d4[1,3])  #nearc value 4 - #nearc value 1 /nearc value_ne 4 - #nearc value 1
# 0.21 of increase nearc 

(d4[4,5] -d4[3,3])/ (d4[4,3] -d4[3,3])  #palearc value 4 - #palearc value 1 /palearc value_ne 4 - #palearc value 1
# 0.68 of increase palearc

(d4[6,5] -d4[5,3])/ (d4[6,3] -d4[5,3])  #tundra value 4 - #tundra value 1 /tundra value_ne 4 - #tundra value 1
# 0.5 of decrease tundra

