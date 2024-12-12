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

#### data ####
setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM")

path <- "C:/Users/roschw001/Documents/migration values"
traits <- read.csv("C:/Users/roschw001/Documents/migration values/migration_meters_all.csv")
h0 <- read.csv("hypotheses_results/H0_new/results_sigma1new.csv", sep=";", dec=",")

h0t <- h0 %>%
  left_join(traits, by = "species")
h0t$gform[134] <- "Shrub"
h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"

result <- as.data.frame(matrix(nrow = 1174, ncol=3))
result[,1] <- h0t$species
result[,2] <- h0t$mode
result[,3] <- h0t$gform
result[,4] <- h0t$class
result[,5] <- h0t$pres

colnames(result) <- c("species", "mode", "gform","class","cells")

result <- result %>%
  mutate(mode = case_when(
    mode == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(mode)  # Keep the original value as a string if none of the above conditions match
  ))

# result <- result %>%
#   mutate(gform = case_when(
#     gform == "aquatic" ~ "Aquatic",
#     TRUE ~ as.character(mode)  # Keep the original value as a string if none of the above conditions match
#   ))

###### mode ####

result1 <- result %>%
  count(mode, name = "total") %>%
  mutate(percentage = (total / sum(total)) * 100) %>%
  arrange(desc(total)) 

##### class ####

result2 <- result %>%
  count(class, name = "total") %>%
  mutate(percentage = (total / sum(total)) * 100) %>%
  arrange(desc(total)) 

##### gform ####

result3 <- result %>%
  count(gform, name = "total") %>%
  mutate(percentage = (total / sum(total)) * 100) %>%
  arrange(desc(total)) 

result3 <- result3 %>%
  mutate(gform = case_when(
    gform == "Submerged/floating aquatic" ~ "aquatic",
    TRUE ~ as.character(gform)  # Keep the original value as a string if none of the above conditions match
  ))

##### arrange data ####

dat <- as.data.frame(cbind(result$gform, result$class))

colnames(dat) <- c("node", "next_node")

# Count occurrences for each flow
flow_data <- as.data.frame(table(dat))

# Rename columns for clarity
colnames(flow_data) <- c("node", "next_node", "count")


flow_data <- flow_data  %>%
  mutate(node = case_when(
    node == "Submerged/floating aquatic" ~ "Aquatic",
    TRUE ~ as.character(node)  # Keep the original value as a string if none of the above conditions match
  ))


flow_data <- flow_data  %>%
  mutate(next_node = case_when(
    next_node == "1" ~ "1: 1 m",
    next_node == "2" ~ "2: 5 m",
    next_node == "3" ~ "3: 15 m",
    next_node == "4" ~ "4: 150 m",
    next_node == "5" ~ "5: 500 m",
    next_node == "6" ~ "6: 1500 m",
    TRUE ~ as.character(node)  # Keep the original value as a string if none of the above conditions match
  ))

# Create the Sankey plot
# Load libraries
library(ggplot2)
library(ggalluvial)
library(viridis) 

#### sankey ####
sankey_plot <- ggplot(flow_data,
                      aes(axis1 = node,
                          axis2 = next_node,
                          y = count)) +
  geom_alluvium(aes(fill = node)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_fill_viridis_d(option = "D") +  # Use a color-blind-friendly 
  #scale_fill_brewer(type = "qual") +
  theme_minimal() +
  labs(
    x = "",
    y = "Number of species",
    title = "B)")+
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),    # Remove x-axis text
    axis.ticks.x = element_blank(),   # Remove x-axis ticks
    axis.line.x = element_blank(),      # Remove x-axis line
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
  legend.position = "none",
  options(repr.plot.width = 5, repr.plot.height = 7) )           # Remove legend

# Display the plot
print(sankey_plot)

ggsave("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/plots/f1b.png", sankey_plot, units = "cm", width = 14, height = 14)


