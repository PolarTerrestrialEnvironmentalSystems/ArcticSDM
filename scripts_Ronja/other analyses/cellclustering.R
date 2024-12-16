#### cellclustering ####

#source louvain_analysis
#https://stackoverflow.com/questions/49834827/louvain-community-detection-in-r-using-igraph-format-of-edges-and-vertices

#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)

library(psych)
library(igraph)

#### clustering ####

data <- "/home/pd/roschw001/R/data"

load("/bioing/user/roschw/sdm_occurance_metadata_25km_dispersal_99.rda") #grid
load(glue::glue("{data}/sigma1Array.rda")) #array restricted dispersal

#time series
for(y in 1:4){
dat <- spArray[,,y,3] #year

dat[is.na(dat)] <- 0
cormatrix <- cor (dat)

distancematrix <- cor2dist(cormatrix)
    
DM2 <- as.matrix(distancematrix)

DM2[cormatrix < 0.3] = 0
DM2[is.na(DM2)] <- 0

G2a <- graph_from_adjacency_matrix(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)
cl2 <- cluster_louvain(G2a) 

membercl <- data.frame(group = cl2$membership, label = 1:38657)

write.csv(membercl, glue::glue("/bioing/user/roschw/membercl_cells_y{y}s3_sigma1.csv"), row.names = FALSE)


#### plotting ####
# Step 1: Identify groups with only one label
single_label_groups <- names(which(table(membercl$group) == 1))

# Step 2: Update the group values to 0 for these labels
membercl$group[membercl$group %in% as.numeric(single_label_groups)] <- 0

#load grid
data <- "/bioing/user/roschw/"
load(glue::glue("{data}/sdm_occurance_metadata_25km_dispersal_99.rda")) #grid
grid <- as.data.frame(spArrayMeta[[2]])

com <- grid %>% mutate(sp = membercl$group) %>% #community
  st_as_sf()

library(ggplot2)

g <- ggplot(com) +
  geom_sf(aes(color = as.factor(sp)), size = 0.001) +
  theme_minimal() +
  labs(title = glue::glue("y{y}s3 sigma1 dispersal"), color = "cell communities")


ggsave(glue::glue("/bioing/user/roschw/cellcommunities_y{y}s3_sigma1.png"), plot = g)

}
