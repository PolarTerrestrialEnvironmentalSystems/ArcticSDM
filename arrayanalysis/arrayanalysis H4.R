#arrayanalysis H2

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
    
    


#communities
subset(communities[[1]][[1]], communities[[1]][[1]][1]==1)


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

cl <- cluster_louvain(G2)
plot(G2, vertex.color=rainbow(4, alpha=0.6)[cl$membership])


cl$names
cl$membership
membercl <- data.frame(group = cl$membership, label = cl$names)

communities[[y]][[s]] <- membercl


}
}

#### community tab list ####
#4y scen list
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


