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

#input
#setwd("//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays")
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays"

load(glue::glue("{data}/spTable.rda"))
load(glue::glue("{data}/binaryArray_25km.rda"))

#binarArray[species, cells, year, scenario]


#all species y1s1
sppres <- as.data.frame(binarArray[,,1,1])
dat <- as.data.frame(t(sppres))

site <- seq(1, nrow(dat),1)

colnames(dat) <- spTable$species

cormatrix <- as.data.frame(cor(dat))

data <- cbind(site, dat)

d <- as.data.frame(t(data))

#### H4 net y1s1 ####
library(psych)
library(igraph)

sppres <- as.data.frame(binarArray[,,1,1])
dat <- as.data.frame(t(sppres))
colnames(dat) <- spTable$species
    
cormatrix <- as.data.frame(cor(dat))
distancematrix <- cor2dist(cormatrix)
    
DM2 <- as.matrix(distancematrix)
## Zero out connections where there is low correlation
DM2[cormatrix < 0.3] = 0
    
G2 <- graph.adjacency(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)
    
vcount(G2)
ecount(G2)
    
library(igraph)

G = G2 #graph

# 1. Apply clustering algorithm A on G nP times, yielding nP partitions
nP <- 10 # number of partitions to generate
partitions <- list()
for (i in 1:nP) {
  partitions[[i]] <- cluster_fast_greedy(G) # replace with desired clustering algorithm
}


# 2. Compute the consensus matrix D
n <- vcount(G)
D <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    same_community <- logical(nP)
    for (p in 1:nP) {
      same_community[p] <- length(intersect(partitions[[p]]$names[which(partitions[[p]]$membership == partitions[[p]]$membership[i])], partitions[[p]]$names[which(partitions[[p]]$membership == partitions[[p]]$membership[j])])) > 0
    }
    D[i,j] <- sum(same_community) / nP
  }
}


# 3. Set entries below threshold τ to zero
tau <- 0.5 
D[D < tau] <- 0

# 4. Apply clustering algorithm A on D nP times, yielding nP partitions
partitions_D <- list()
for (i in 1:nP) {
  partitions_D[[i]] <- cluster_louvain(graph_from_adjacency_matrix(D, mode = "undirected")) 
}


# 5. Check if partitions are equal
partitions_equal <- TRUE
for (i in 2:nP) {
  if (!all(partitions_D[[1]]$membership == partitions_D[[i]]$membership)) {
    partitions_equal <- FALSE
    break
  }
}


if (!partitions_equal) {
  # Repeat from step 2
  D <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (p in 1:nP) {
        if (partitions_D[[p]][i] == partitions_D[[p]][j]) {
          D[i,j] <- D[i,j] + 1
        }
      }
      D[i,j] <- D[i,j] / nP
    }
  }
  D[D < tau] <- 0
  partitions_D <- list()
  for (i in 1:nP) {
    partitions_D[[i]] <- cluster_fast_greedy(graph_from_adjacency_matrix(D))
  }
  
  # Check again
  partitions_equal <- TRUE
  for (i in 2:nP) {
    if (!all(partitions_D[[1]] == partitions_D[[i]])) {
      partitions_equal <- FALSE
      break
    }
  }
}

if (partitions_equal) {
  print("Clustering converged")
  print(partitions_D[[1]])
} else {
  print("Clustering did not converge")
}

membercl <- data.frame(group = partitions_D[[1]]$membership, label = spTable$species)
membercl #communities for y1s1
    
  
#### H5 net loop ####
library(psych)
library(igraph)

communities <- list()

for (y in 1:4){
  communities[[y]] <- list()
  for (s in 1:3){  

sppres <- as.data.frame(binarArray[,,y,s])
dat <- as.data.frame(t(sppres))
colnames(dat) <- spTable$species

cormatrix <- as.data.frame(cor(dat))
distancematrix <- cor2dist(cormatrix)

DM2 <- as.matrix(distancematrix)
## Zero out connections where there is low correlation
DM2[cormatrix < 0.3] = 0

G2 <- graph.adjacency(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)

vcount(G2)
ecount(G2)

##### consensus matrix perplexity ####

library(igraph)

G = G2 #graph

# 1. Apply clustering algorithm A on G nP times, yielding nP partitions
nP <- 10 # number of partitions to generate
partitions <- list()
for (i in 1:nP) {
  partitions[[i]] <- cluster_fast_greedy(G) # replace with desired clustering algorithm
}


# 2. Compute the consensus matrix D
n <- vcount(G)
D <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    same_community <- logical(nP)
    for (p in 1:nP) {
      same_community[p] <- length(intersect(partitions[[p]]$names[which(partitions[[p]]$membership == partitions[[p]]$membership[i])], partitions[[p]]$names[which(partitions[[p]]$membership == partitions[[p]]$membership[j])])) > 0
    }
    D[i,j] <- sum(same_community) / nP
  }
}


# 3. Set entries below threshold τ to zero
tau <- 0.5 
D[D < tau] <- 0

# 4. Apply clustering algorithm A on D nP times, yielding nP partitions
partitions_D <- list()
for (i in 1:nP) {
  partitions_D[[i]] <- cluster_louvain(graph_from_adjacency_matrix(D, mode = "undirected")) 
}


# 5. Check if partitions are equal
partitions_equal <- TRUE
for (i in 2:nP) {
  if (!all(partitions_D[[1]]$membership == partitions_D[[i]]$membership)) {
    partitions_equal <- FALSE
    break
  }
}


if (!partitions_equal) {
  # Repeat from step 2
  D <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (p in 1:nP) {
        if (partitions_D[[p]][i] == partitions_D[[p]][j]) {
          D[i,j] <- D[i,j] + 1
        }
      }
      D[i,j] <- D[i,j] / nP
    }
  }
  D[D < tau] <- 0
  partitions_D <- list()
  for (i in 1:nP) {
    partitions_D[[i]] <- cluster_fast_greedy(graph_from_adjacency_matrix(D))
  }
  
  # Check again
  partitions_equal <- TRUE
  for (i in 2:nP) {
    if (!all(partitions_D[[1]] == partitions_D[[i]])) {
      partitions_equal <- FALSE
      break
    }
  }
}

if (partitions_equal) {
  print("Clustering converged")
  print(partitions_D[[1]])
} else {
  print("Clustering did not converge")
}


membercl <- data.frame(group = partitions_D[[1]]$membership, label = spTable$species)
communities[[y]][[s]] <- membercl

  }
}



#communities
communities[[1]][[3]]
communities[[1]][[1]]

#### community tab list ####
#4y scen list

#communities[[year]][[scen]][group]
comlist <- list()

for (s in 1:3){
  d <-  as.data.frame(cbind(communities[[1]][[s]][2], communities[[1]][[s]][1], 
                            communities[[2]][[s]][1], communities[[3]][[s]][1], communities[[4]][[s]][1]))
  colnames(d) <- c("species", "y1", "y2", "y3", "y4")
  
  comlist[[s]] <- d
}

comlist[[1]]

#### scenario tab list ####

scenlist <- list()

for (y in 1:4){
  d <-  as.data.frame(cbind(communities[[y]][[1]][2], communities[[y]][[1]][1], 
                            communities[[y]][[2]][1], communities[[y]][[3]][1]))
  colnames(d) <- c("species", "s1", "s2", "s3")
  
  scenlist[[y]] <- d
}

scenlist[[1]]


#### jaccard ####

#Calculate Jaccard similarity
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(unique(c(a, b)))
  jac <- intersection / union
}

j <- as.data.frame(matrix(ncol=4, nrow=3))
colnames(j) <- c(1,2,3,4)
rownames(j) <- c(1,2,3)


jlist <- list()

for (c in 1:4){ #c is the community number
  
a <- subset(communities[[1]][[1]], communities[[1]][[1]][1]==c)

for (y in 1:4){
  for (s in 1:3){

    b <- subset(communities[[y]][[s]], communities[[y]][[s]][1]==c)
    
j[s,y] <- jaccard(a,b)

  }
}
  jlist[[c]] <- j
}


jlist[[4]]




