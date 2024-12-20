#1 with table 12 rows
#2 with table 4 rows, per scenario
#3 with species list and tab, 12 rows
#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(ggplot2)

#### general ####
setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM")
path <- "C:/Users/roschw001/Documents/migration values"

traits <- read.csv("C:/Users/roschw001/Documents/migration values/migration_meters_all.csv")
#h0 <- read.csv("hypotheses_results/H0/H0_percenttab_unres.csv", sep=";", dec=",")


#### p1 ####
h0 <- read.csv("hypotheses_results/H0_new/results_sigma1new.csv", sep=";")
h0t <- h0 %>%
  left_join(traits, by = "species")

h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"
h0t$gform[134] <- "Shrub"

##### mode ####

result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$mode
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "mode", "new_cells")

result <- result %>%
  mutate(mode = case_when(
    mode == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(mode)  # Keep the original value as a string if none of the above conditions match
  ))

result1 <-  result

result1a <- subset(result1, result1$new_cells > 0)
sum(result1a$mode =="Myrmecochory")
sum(result1$mode =="Myrmecochory")

result1a <- subset(result1, result1$new_cells > 0)
sum(result1a$mode =="local")
sum(result1$mode =="local")

result1 <- subset(result1, result1$mode !="local")

# Create a boxplot
p1 <- ggplot(result1, aes(y = reorder(mode, -new_cells, FUN = median), x = new_cells)) +
  geom_boxplot(fill = "#0072B2") +
  labs(#title = expression(paste(Delta, " 2100 constrained")),  # Update labels    
       y = "dispersal mode",
       x = "new cells") +
  theme_minimal()+
  ggtitle("B)")+
  theme(margin(t = 10, r = 10, b = 60, l = 10))+
  theme(   plot.title = element_text(size = 16),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

print(p1)
##### class ####

result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$class
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "class", "new_cells")

result <- result %>%
  mutate(class = case_when(
    class == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(class)  # Keep the original value as a string if none of the above conditions match
  ))

#result1m <- subset(result, result$new_cells > 0)
result1m <-  result

result1a <- subset(result1m, result1m$new_cells > 0)
sum(result1a$class == "2")
sum(result1m$class == "2")

result1a <- subset(result1m, result1$new_cells > 0)
sum(result1a$mode =="local")
sum(result1$mode =="local")

result1m$classn <- as.numeric(result1m$class)
result1m <-  subset(result1m, result1m$classn > 2)

p1m <- ggplot(result1m, aes(y = reorder(class, -new_cells, FUN = median), x = new_cells)) +
  geom_boxplot(fill = "#D95F02") +
  labs(title = "",
       y = "dispersal class",
       x = "new cells") +
  theme_minimal()+
  ggtitle("C)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 18))
      
print(p1m)
##### gform ####

result<- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$gform
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "gform", "new_cells")
unique(result$gform)

result <- result %>%
  mutate(gform = case_when(
    gform == "Submerged/floating aquatic" ~ "aquatic",
    TRUE ~ as.character(gform)  # Keep the original value as a string if none of the above conditions match
  ))

result1g <- result #subset(result, result$new_cells > 0)
# Load ggplot2 library
library(ggplot2)

# Create a boxplot
p1g <- ggplot(result1g, aes(y = reorder(gform, -new_cells, FUN = median), x = new_cells)) +
  geom_boxplot(fill = "#009E73") +
  labs(title = "",
       y = "growth form",
       x = "new cells") +
  theme_minimal()+
  ggtitle("A)")+
  theme(plot.title = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))
print(p1g)

#### p2 ####
h0 <- read.csv("hypotheses_results/H0_new/results_sig5.csv", sep=";")

h0t <- h0 %>%
  left_join(traits, by = "species")

h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"
h0t$gform[134] <- "Shrub"

##### mode ####

result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$mode
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "mode", "new_cells")

result <- result %>%
  mutate(mode = case_when(
    mode == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(mode)  # Keep the original value as a string if none of the above conditions match
  ))

result2 <-  result #subset(result, result$new_cells > 0)

p2 <- ggplot(result2, aes(x = reorder(mode, -new_cells, FUN = median), y = new_cells)) +
  geom_boxplot() +
  labs(title = "sig5",
       x = "dispersal mode",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



result_mode <- left_join(result1, result2, by = "species")
result_mode$ratio <- result_mode$new_cells.x/result_mode$new_cells.y

pmode <- ggplot(result_mode, aes(x = reorder(mode.x, -ratio, FUN = median), y = ratio)) +
  geom_boxplot() +
  labs(title = "sig1/sig5",
       x = "dispersal mode",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot(pmode)

##### class ####

result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$class
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "class", "new_cells")

result <- result %>%
  mutate(class = case_when(
    class == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(class)  # Keep the original value as a string if none of the above conditions match
  ))

result2m <-  result #subset(result, result$new_cells > 0)

p2m <- ggplot(result2m, aes(x = reorder(class, -new_cells, FUN = median), y = new_cells)) +
  geom_boxplot() +
  labs(title = "sig5",
       x = "dispersal class",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

result_class <- left_join(result1m, result2m, by = "species")
result_class$ratio <- result_class$new_cells.x/result_class$new_cells.y

pclass <- ggplot(result_class, aes(x = reorder(class.x, -ratio, FUN = median), y = ratio)) +
  geom_boxplot() +
  labs(title = "sig1/sig5",
       x = "dispersal class",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot(pclass)

##### gform ####
result<- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$gform
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "gform", "new_cells")
unique(result$gform)

result <- result %>%
  mutate(gform = case_when(
    gform == "Submerged/floating aquatic" ~ "aquatic",
    TRUE ~ as.character(gform)  # Keep the original value as a string if none of the above conditions match
  ))

result2g <-  result #subset(result, result$new_cells > 0)
# Load ggplot2 library
library(ggplot2)

# Create a boxplot
p2g <- ggplot(result2g, aes(x = reorder(gform, -new_cells, FUN = median), y = new_cells)) +
  geom_boxplot(fill = "#009E73") +
  labs(title = "sig5",
       x = "Growth form",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p2g)


result1g$sig5 <- result2g$new_cells 
result1g$ratio <- result1g$new_cells/result1g$sig5

# Create a boxplot
p12g <- ggplot(result1g, aes(x = reorder(gform, -ratio, FUN = median), y = ratio)) +
  geom_boxplot(fill = "#009E73") +
  labs(title = "sig1/sig5",
       x = "Growth form",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p12g)

#### p3 ####
h0 <- read.csv("hypotheses_results/H0_new/results_unres.csv", sep=";")
h0 <- h0[-341,]
h0t <- h0 %>%
  left_join(traits, by = "species")

h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"
h0t$gform[134] <- "Shrub"

##### mode ####
result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$mode
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "mode", "new_cells")

result <- result %>%
  mutate(mode = case_when(
    mode == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(mode)  # Keep the original value as a string if none of the above conditions match
  ))

result3 <-  result #subset(result, result$new_cells > 0)
 
p3 <- ggplot(result3, aes(x = reorder(mode, -new_cells, FUN = median), y = new_cells)) +
  geom_boxplot() +
  labs(title = "unres",
       x = "dispersal mode",
       y = "new cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

##### class ####

result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$class
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "class", "new_cells")

result <- result %>%
  mutate(class = case_when(
    class == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(class)  # Keep the original value as a string if none of the above conditions match
  ))

result3m <-  result #subset(result, result$new_cells > 0)

p3m <- ggplot(result3m, aes(x = reorder(class, -new_cells, FUN = median), y = new_cells)) +
  geom_boxplot() +
  labs(title = "unres",
       x = "dispersal class",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



##### gform ####
result<- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$gform
result[,3] <- h0t$y4s3

colnames(result) <- c("species", "gform", "new_cells")
unique(result$gform)

result <- result %>%
  mutate(gform = case_when(
    gform == "Submerged/floating aquatic" ~ "aquatic",
    TRUE ~ as.character(gform)  # Keep the original value as a string if none of the above conditions match
  ))

result3g <-  result #subset(result, result$new_cells > 0)
# Load ggplot2 library
library(ggplot2)

# Create a boxplot
p3g <- ggplot(result3g, aes(x = reorder(gform, -new_cells, FUN = median), y = new_cells)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "unres",
       x = "Growth form",
       y = "new cells") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p3)

#### combined ####
# library(cowplot)
# combined_plot <- plot_grid(p1, p3, p2, ncol = 3)
# print(combined_plot)
# 
# combined_plot_g <- plot_grid(p1g, p3g, p2g, ncol = 3)
# print(combined_plot_g)
# 
# combined_plot_m <- plot_grid(p1m, p3m, p2m, ncol = 3)
# print(combined_plot_m)
##############################################


# combined_plot_both <- plot_grid(p1m, p1, ncol = 2)
# print(combined_plot_both)
# 
# combined_plot_both <- plot_grid(
#   p1m + theme(plot.title = element_blank()),
#   p1 + theme(plot.title = element_blank()),
#   p1g + theme(plot.title = element_blank()),
#   ncol = 3
# )
# 
# print(combined_plot_both)


library(cowplot)

newcells <- plot_grid(p1g, p1, p1m, ncol=3)
print(newcells)


combined_plot_both <- plot_grid(
  p1m + theme(plot.title = element_blank()),
  p1 + theme(plot.title = element_blank()),
  p1g + theme(plot.title = element_blank()),
  ncol = 2
)


print(combined_plot_both)


combined <- plot_grid(combined_plot_both, combined_plot_g, nrow = 2)
print(combined)


################################################

# library(cowplot)
# 
# # Suppress y-axes for each plot
# p1m0 <- p1m +  theme(plot.title = element_blank())
# 
# # + #theme(axis.title.y = element_blank())
#                    #axis.text.y = element_blank())
#                    #axis.ticks.y = element_blank())
# 
# p10 <- p1 + theme(axis.title.y = element_blank(), 
#                  axis.text.y = element_blank(), 
#                  axis.ticks.y = element_blank(),
#                  plot.title = element_blank())
# 
# p1g0 <- p1g + theme(axis.title.y = element_blank(), 
#                    axis.text.y = element_blank(), 
#                    axis.ticks.y = element_blank(),
#                    plot.title = element_blank())
# 
# 
# # Add a new y-axis label
# combined_plot_both <- ggdraw() +
#   draw_label("New Y-Axis Label", angle = 90, vjust = 0.5, x = -0.05) +   # Adjust x for positioning
#   draw_plot(combined_plot_both)
# 
# 
# combined_plot_both <- plot_grid(p1m0, p10, p1g0, ncol = 3)
# print(combined_plot_both)
# 
# # Create a dummy plot for the y-axis
# y_axis_plot <- ggplot() +
#   scale_y_continuous(breaks = seq(0, 6000, by = 1000)) +  # Set y-ticks from 0 to 6000
#   labs(y = "New Y-Axis Label") +                         # Add y-axis label
#   theme_minimal() +
#   theme(axis.title.y = element_text(angle = 90, vjust = 0.5), # Rotate label
#         axis.text.y = element_text(size = 10),          # Adjust text size
#         axis.ticks.y = element_line(size = 0.5),        # Adjust tick size
#         panel.grid.major.y = element_line(color = "grey"), # Optional: Add grid lines
#         panel.grid.minor.y = element_blank(),            # Remove minor grid lines
#         panel.grid.major.x = element_blank(),            # Remove x grid lines
#         axis.title.x = element_blank(),                   # Remove x title
#         axis.text.x = element_blank(),                    # Remove x text
#         axis.ticks.x = element_blank())                   # Remove x ticks
# 
# print(y_axis_plot)
# 
# # Combine with the dummy y-axis plot on the left
# final_combined_plot <- plot_grid(
#   y_axis_plot,
#   combined_plot_both,
#   ncol = 2,
#   rel_widths = c(0.1, 1)    # Adjust width ratio as needed
# )
# 
# # Display the final combined plot
# print(final_combined_plot)
# 
# 
# unique(result1$mode)
# sum(result1$mode =="Myrmecochory")
# sum(result1$mode =="local")
# subset(result1, result1$mode =="Myrmecochory")
# sum(result1$mode =="local")
#        