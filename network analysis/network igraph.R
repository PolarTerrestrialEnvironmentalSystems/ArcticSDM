#network


install.packages(c("igraph", "tidygraph", "ggraph", "quanteda"))
library(igraph)
library(tidyr)
library(tidygraph)
library(ggraph)
library(quanteda)

#input

library(stars)
load("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/occTab_all.rda")

species <- occTab_all$species
coordsxy <- as.data.frame(st_coordinates(occTab_all, dims = c("x", "y")))

coordsxy$x1 = round(coordsxy$X, -4)
coordsxy$y1 = round(coordsxy$Y, -4)
nettab = cbind(coordsxy[,3:4], species)


# Unite 'col1' and 'col2' into a single column 'united_col'
nettab <- unite(nettab, col = site, x1, y1, sep = "_")


# Example species observations data
species_obs <- data.frame(
  species = nettab$species,
  location = nettab$site
)

write.csv(species_obs, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/speices_obs.csv")
getwd()
so <- species_obs[1:10000,]
net <- graph_from_data_frame(so, directed = FALSE)

write.csv(species_obs, "C:/Users/roschw001/Documents/R/SDM/SDM/data nets.csv")
# Explore the network
V(net)
E(net)

# Visualize the network
plot(net)

net2 <- simplify(net)

plot(net2, vertex.label = NA)

plot(net, edge.color="gray", vertex.color="lightblue", vertex.label = NA) 

#### heat map ####

netm <- as_adjacency_matrix(net2)

colnames(netm) <- V(net)$

rownames(netm) <- V(net)$'List of 4'$names

V(net)$'List of 4'$'List of 20'

palf <- colorRampPalette(c("gold", "dark orange")) 

heatmap(netm, Rowv = NA, Colv = NA, col = palf(100), 
        
        scale="none", margins=c(10,10) )


plot(closeness(net, mode="all", weights=NA))

plot(cocitation(net))

net.sym <- as.undirected(net, mode= "collapse",
                         
                         edge.attr.comb=list(weight="sum", "ignore"))

#### Find cliques ####(complete subgraphs of an undirected graph)

clique <- cliques(net.sym) # list of cliques       

sapply(cliques(net.sym), length) # clique sizes

largest_cliques(net.sym) # cliques with max number of nodes

#### community detection ####

ceb <- cluster_edge_betweenness(net) 
ceb2 <- cluster_walktrap(net)
ceb3 <- cluster_fast_greedy(net)
ceb4 <- cluster_leading_eigen(net)

dendPlot(ceb, mode="hclust")

plot(ceb, net, vertex.label=NA)
plot(ceb2, net, vertex.label=NA)
plot(ceb3, net, vertex.label=NA)
plot(ceb4, net, vertex.label=NA)
ceb


edge.betweeness.community(net)

ceb$bridges
class(ceb)

## [1] "communities"

length(ceb)     # number of communities

## [1] 5

membership(ceb) # community membership for each node

## s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 s11 s12 s13 s14 s15 s16 s17 

##   1   2   3   4   1   4   3   3   5   5   4   4   4   4   1   4   4

modularity(ceb) # how modular the graph partitioning is

## [1] 0.292476

crossing(ceb, net)   # boolean vector: TRUE for edges across communities



#### plot types ####
fg <- make_full_graph(net)

plot(fg, vertex.size=10, vertex.label=NA)


st <- make_star(40)

plot(st, vertex.size=10, vertex.label=NA) 

#Tree graph

tr <- make_tree(40, children = 3, mode = "undirected")

plot(tr, vertex.size=10, vertex.label=NA) 

#Erdos-Renyi random graph model
# (‘n’ is number of nodes, ‘m’ is the number of edges).

er <- sample_gnm(n=100, m=40) 

plot(er, vertex.size=6, vertex.label=NA)  

#Watts-Strogatz small-world model
sw <- sample_smallworld(dim=2, size=10, nei=1, p=0.1)

plot(sw, vertex.size=6, vertex.label=NA, layout=layout_in_circle)

# Calculate degree centrality
degree_centrality <- igraph::degree(net)

# Calculate betweenness centrality
betweenness_centrality <- igraph::betweenness(net)

# Calculate clustering coefficient
clustering_coef <- igraph::transitivity(net)

# Create a ggraph object
g <- ggraph(net, layout="fr")

ggraph(g) + 
  geom_edge_link0(aes(color=as.factor(color)), width=0.6, alpha=0.35) + 
  geom_node_point(aes(color=as.factor(color)), size=3, alpha=0.75) + 
  theme_graph(base_family = 'Helvetica')

ggraph(g, layout='fr') + 
  geom_edge_link0(aes(filter=color!=9999 ,color=as.factor(color)), width=0.6, alpha=0.35) + 
  geom_edge_link0(aes(filter=color==9999), color='grey', width=0.5, alpha=0.25) + 
  geom_node_point(aes(color=as.factor(color)), size=3, alpha=0.75) + 
  theme_graph(base_family = 'Helvetica')

# Visualize the network
g + geom_edge_link() + geom_node_point()

summary(net)


##walktrap
library(igraph)
all_wt<- walktrap.community(net, steps=6,modularity=TRUE)
all_wt_memb <- community.to.membership(all, all_wt$merges, steps=which.max(all_wt$modularity)-1)


colbar <- rainbow(20)
col_wt<- colbar[all_wt_memb$membership+1]

l <- layout.fruchterman.reingold(all, niter=100)
plot(all, layout=l, vertex.size=3, vertex.color=col_wt, vertex.label=NA,edge.arrow.size=0.01,
     main="Walktrap Method")

