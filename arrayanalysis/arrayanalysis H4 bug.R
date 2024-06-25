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
    
    cl <- cluster_louvain(G2)
    plot(G2, vertex.color=rainbow(4, alpha=0.6)[cl$membership])
    
    
    cl$names
    cl$membership
    membercl <- data.frame(group = cl$membership, label = cl$names)
    
    membercl #communities for y1s1
    
    
#sppresdif <- sppres3-sppres1

#### H5 net loop ####
library(psych)
library(igraph)

communities <- list()

for (y in 1:4){
  communities[[y]] <- list()
  for (s in 1:3){  
 
y = 1
s = 1

sppres <- as.data.frame(binarArray[,,y,s])
dat <- as.data.frame(t(sppres3))
colnames(dat) <- spTable$species

cormatrix <- as.data.frame(cor(dat))
#distancematrixdif <- distancematrix3 - distancematrix1
distancematrix <- cor2dist(cormatrix)

DM2 <- as.matrix(distancematrix)
## Zero out connections where there is low correlation
DM2[cormatrix < 0.3] = 0

#Gdif <- as.data.frame(G23 - G21)
#DMdif
G2 <- graph.adjacency(DM2, mode = "undirected", weighted = TRUE, diag = TRUE)

vcount(G2)
ecount(G2)

plot(G2, vertex.color=rainbow(4, alpha=0.6)[cl$membership])

cl$names
cl1$membership
membercl <- data.frame(group = cl$membership, label = cl$names)
membercl
communities[[y]][[s]] <- membercl


  }
}

#### consensus matrix perplexity ####

library(igraph)

G = G2 #graph

# 1. Apply clustering algorithm A on G nP times, yielding nP partitions
nP <- 16 # number of partitions to generate
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

partitions_D[[16]] #same for 1-16


#### iris geek ####

#Using Euclidean distance as a similarity measure
similarity_matrix <- exp(-dist(dat)^2 / (2 * 1^2))

#Compute Eigenvalues and Eigenvectors
eigen_result <- eigen(similarity_matrix)
eigenvalues <- eigen_result$values
eigenvectors <- eigen_result$vectors

#Choose the First k Eigenvectors
k <- 3 
selected_eigenvectors <- eigenvectors[, 1:k]

#Apply K-Means Clustering
cluster_assignments <- kmeans(selected_eigenvectors, centers = k)$cluster

# Add species information to the clustering results
iris$Cluster <- factor(cluster_assignments)
iris$Species <- as.character(iris$Species)


library(ggplot2)

# Visualizing the clusters with species names
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Cluster, label = Species)) +
  geom_point() +
  geom_text(check_overlap = TRUE, vjust = 1.5) +
  labs(title = "Spectral Clustering of Iris Dataset",
       x = "Sepal Length", y = "Sepal Width")


#### spectral clustering perplexity ####

library(igraph)
library(kernlab)

set.seed(123) # Set a seed for reproducibility

# For a graph
g <- G2
sim_matrix <- as_adjacency_matrix(g, attr="weight", sparse=FALSE)

k <- 3  # Number of communities/clusters

# library(RSpectra)
# dim(sim_matrix)
# 
# ev <- arpack(sim_matrix, k)
# ev <- eigs_sym(sim_matrix, k)
# communities <- ev$vectors[, 1]

sc <- specc(sim_matrix, centers=k, nystrom.red=FALSE, initseed=123)
communities <- sc@.Data

communities <- sort(communities)
communities

# Set vertex colors based on community assignments
V(g)$color <- rainbow(k)[communities]

# Plot the graph with colored vertices
plot(g, vertex.label=g$name)


library(FactoMineR)

# Perform PCA
res.pca <- PCA(dat, scale.unit = TRUE, graph = FALSE)

# Extract principal component scores
pca_scores <- res.pca$ind$coord

# Perform k-means clustering on PCA scores
k <- 3  # Number of clusters
km_res <- kmeans(pca_scores, centers = k)
km_res$cluster
#### Spectral Clustering with k-means geeks ####

# Compute the similarity matrix
similarity_matrix <- exp(-dist(dat)^2)

# Perform spectral decomposition
eigen_result <- eigen(similarity_matrix)

# Extract the top-k eigenvectors
k_eigenvectors <- eigen_result$vectors[, 1:k]

# Perform k-means clustering on the eigenvectors
cluster_assignments <- kmeans(k_eigenvectors, centers = k)$cluster

# Visualize the clusters
plot(data, col = cluster_assignments, pch = 19, 
     main = "Spectral Clustering with k-means")



#### cvx cluster ####

library(cvxclustr)
create_clustering_problem(dat)
cxv <- cvxclust_ama(dat, nu=1/ncol(dat))
cvxclust_admm()
cl <- cluster_louvain(G2)
cl <- cluster_edge_betweenness(G2)
cl <- cluster_fast_greedy(G2)


#### 


#communities
subset(communities[[1]][[3]], communities[[1]][[3]][1]==1)


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

#communities[[y]][[s]]
y1s1<- subset(communities[[1]][[1]], communities[[1]][[1]][1]==1)
y3s1 <- subset(communities[[4]][[1]], communities[[4]][[1]][1]==1)

communities[[1]][[1]][1]==2
#Calculate Jaccard similarity
# library(stringdist)
# jaccard_sim <- jacc_dist(y1s1, y3s1)
# vegdist(x, method="jaccard", binary = T)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(unique(c(a, b)))
  jac <- intersection / union
}

j <- as.data.frame(matrix(ncol=4, nrow=3))
colnames(j) <- c(1,2,3,4)
rownames(j) <- c(1,2,3)


jlist <- list()

for (c in 1:4){
  
a <- subset(communities[[1]][[1]], communities[[1]][[1]][1]==c)

for (y in 1:4){
  for (s in 1:3){

    b <- subset(communities[[y]][[s]], communities[[y]][[s]][1]==c)
    
j[s,y] <- jaccard(a,b)

  }
}
  jlist[[c]] <- j
}


jlist[[2]]
#communities[[y]][[s]]
communities[[1]][[3]]


