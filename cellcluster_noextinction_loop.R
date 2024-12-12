#arrayanalysis H2

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

#data <- "/bioing/R/data"
#data <- "/home/pd/roschw001/R/data"

load("//smb.isipd.dmawi.de/projects/bioing/user/roschw/sdm_occurance_metadata_25km_dispersal_99.rda") #grid
load(glue::glue("//smb.isipd.dmawi.de/projects/bioing/user/roschw/sigma1Array.rda")) #array restricted dispersal

y1 <- as.data.frame(spArray[,,1,3])
y2 <- as.data.frame(spArray[,,2,3])
y3 <- as.data.frame(spArray[,,3,3])
y4 <- as.data.frame(spArray[,,4,3])


# Create a function to count the number of new 1 cells for each grid cell

make_comparison<- function(present, future) {
  
  comparison = present + future
  comparison[comparison == 2] <- 1
  return(comparison)
}

present = y1

#### sig1 ####

future <- y4
y1 <- spArray[,,1,3]
y2ne <- make_comparison(spArray[,,1,3], spArray[,,2,3])
y3ne <- make_comparison(spArray[,,1,3], spArray[,,3,3])
y4ne <- make_comparison(spArray[,,1,3], spArray[,,4,3])

# Combine the slices into mixArray
neArray <- abind(y1, y2ne, y3ne, y4ne, along = 3)
a <- as.data.frame(neArray[,,1])

setwd("//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM")
save(neArray, file="neArray.rda")

summary(y4ne)

pres <- colSums(present)
summary(pres)

fut <- colSums(y4)
summary(fut)

both <- colSums(y4ne)
summary(both)

test <- y4ne - y1 #calculate for new cells 
tes <- colSums(test)
summary(tes)


#test
count_new_ones <- function(present, future) {
  new_ones <- future - present
  new_ones[new_ones < 0] <- 0
  return(new_ones)
}


present <- as.data.frame(spArray[, , 1, 1])
future <- as.data.frame(spArray[, , 4, 3])


new_array <- as.data.frame(count_new_ones(present, future))

n1 <- colSums(new_array)

summary(n1)
summary(tes)
# IS THE SAME, calculation CORRECT

sum(n1)
sum(both)
sum(n1)/sum(both)

### only 5.5 % new colonization of cells, not visible
### in the classification, one ore some more species per cells does not make a difference
### need many new species to change a cell's affilation to a cluster

#### clustering loop ####


for(y in 2:4){
  
  future <- spArray[,,y,3] #year
  dat <- as.data.frame(make_comparison(present, future))
  
  dat[is.na(dat)] <- 0
  cormatrix <- cor (dat)
  
  #save(cormatrix, file = "cormatrix_cells_y4s3.rda")
  
  distancematrix <- cor2dist(cormatrix)
  
  DM2 <- as.matrix(distancematrix)
  
  DM2[cormatrix < 0.3] = 0
  DM2[is.na(DM2)] <- 0
  
  G2a <- graph_from_adjacency_matrix(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)
  cl2 <- cluster_louvain(G2a) 
  
  membercl <- data.frame(group = cl2$membership, label = 1:38657)
  
  write.csv(membercl, glue::glue("/bioing/user/roschw/membercl_cells_y{y}s3_ne.csv"), row.names = FALSE)
  
  #save(cl2, file="/bioing/user/roschw/cell_cluster.rda")
  #save(cormatrix, file="/bioing/user/roschw/cormatrix_cells.rda")
  
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
  
  
  ggsave(glue::glue("/bioing/user/roschw/cellcommunities_y{y}s3_ne.png"), plot = g)
  
}

