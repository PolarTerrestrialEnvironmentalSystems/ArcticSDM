# Load required packages
library(factoextra)
library(cluster)
library(stars)
library(tidyr)

load("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/occTab_all.rda")



# Prepare your data (assuming it's in a data frame called 'data')
species <- occTab_all$species
coordsxy <- as.data.frame(st_coordinates(occTab_all, dims = c("x", "y")))

coordsxy$x1 = round(coordsxy$X, -4)
coordsxy$y1 = round(coordsxy$Y, -4)
nettab = cbind(coordsxy[,3:4], species)


# Unite 'col1' and 'col2' into a single column 'united_col'
nettab <- unite(nettab, col = site, x1, y1, sep = "")
nettab$p = 1
# Transform the DataFrame to wide format
netwide <- nettab %>%
  pivot_wider(names_from = species, values_from = p,
              values_fn = list(p = sum)) 



netwide0 <- replace(netwide, is.na(netwide), 0)

library(tidyverse)
#write_csv(netwide0, "C:/Users/roschw001/Documents/R/SDM/SDM/data_nets/netwide0.csv")

#write.csv(netwide0, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/netwide0.csv")

species_obs <- read_csv("C:/Users/roschw001/Documents/R/SDM/SDM/data_nets/species_obs.csv")
speciesobs <- as.data.frame(species_obs[,2:3])
speciesobs$ p = 1
spob <- speciesobs[1:1000,]

netwide0 <- read_csv("C:/Users/roschw001/Documents/R/SDM/SDM/data_nets/netwide0.csv")
netwide2 <- netwide0[1:500,2:50]
netwide2 <- as.data.frame(t(netwide0[1:500,2:50]))

#netwide2 <- replace(netwide1, is.na(netwide1), 0)

# head(netwide2)
# netwide2[1,1]
# x = gsub(",", "", x)
# as.numeric(4830000140000)
# 
# netwide2$site %>% formattable::accounting() %>% as.numeric()
# netwide2$site = as.factor(netwide2$site)
# mean(netwide2$siten)


data_scaled <- scale(netwide2)

# head(data_scaled)
# 
# corsp <- cor(netwide2)
# spob %>%
#   group_by(species, location) %>%
#   cor.test(species, location, method = "spearman")
# 
# head(corsp)

# # Compute Spearman correlation coefficients
# cor_matrix <- cor(netwide2, method = "spearman")
# head(cor_matrix)
# # Compute distance matrix
# dist_matrix <- get_dist(cor_matrix, method = "spearman")

data_sc = as.data.frame(na.omit(data_scaled))

# Perform clustering (e.g., k-means)
clustering_result <- eclust(data_sc, FUNcluster = "kmeans", hc_metric = "spearman")

write_csv(clustering_result, "C:/Users/roschw001/Documents/R/SDM/SDM/data_nets/cluster.csv")

# Visualize clustering results
fviz_cluster(clustering_result, repel = TRUE)

summary(clustering_result)
clustering_result
