#network


#install.packages(c("igraph", "tidygraph", "ggraph", "quanteda"))
library(igraph)
library(tidyr)
library(tidygraph)
library(ggraph)
library(quanteda)
library(stars)

#input

# library(stars)
# load("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/occTab_all.rda")
# 
# species <- occTab_all$species
# coordsxy <- as.data.frame(st_coordinates(occTab_all, dims = c("x", "y")))
# 
# coordsxy$x1 = round(coordsxy$X, -4)
# coordsxy$y1 = round(coordsxy$Y, -4)
# nettab = cbind(coordsxy[,3:4], species)
# 
# 
# # Unite 'col1' and 'col2' into a single column 'united_col'
# nettab <- unite(nettab, col = site, x1, y1, sep = "_")
# 
# 
# # Example species observations data
# species_obs <- data.frame(
#   species = nettab$species,
#   location = nettab$site
# )
# 
# #write.csv(species_obs, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/speices_obs.csv")




netwide0 <- read_csv("C:/Users/roschw001/Documents/R/SDM/SDM/data_nets/netwide0.csv")
netwide2 <- netwide0[1:100,1:1000]

netwide <- netwide2[,2:1000]
rownames(netwide) <- netwide2$site
is.numeric(df_trans)

nw <- box.cox.chord(netwide, bc.exp=0)
df_trans <- nw
############################################################################
# 2 - Co-occurence and community building -> TCorr.test and iGraph
############################################################################

# Do correlation of occurence matrix - very time consuming
corr_df <- corr.test(nw,       # a matrix or dataframe
                     use = "pairwise",  # pairwise" is the default value and will do pairwise deletion of cases. use="complete" will select just complete cases
                     method="spearman", # default is Pearson. spearman is slower especially for large datasets
                     adjust="holm",     # What adjustment for multiple tests should be used? ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
                     alpha=.05,         # alpha level of confidence intervals
                     ci=FALSE,          # default is TRUE but set false. For very large matrices (>200x200), noticeable speed improvement if confidence interval are not found
                     minlength=5,       # What is the minimum length for abbreviations. Defaults to 5.
                     normal=F)          # normal = FALSE -> we do not have normal distribution. But it is much slower to process


############################################################################
# 3 - Prepare for network analysis
############################################################################

# Make p-values as long
corr_p      <- as.data.frame(corr_df$p)
corr_p$sp1  <- rownames(corr_p)
corr_p_long <- gather(corr_p, sp2, prob, all_of(colnames(as.data.frame(corr_df$p)))) 

# Make uniq_ids
corr_p_long$uniqid <- paste(corr_p_long$sp1, corr_p_long$sp2, sep = "_")
corr_p_long        <- corr_p_long %>% select(uniqid, prob)

# Re-format correlation matrix and merge p-value
corr_r      <- as.data.frame(corr_df$r)
corr_r$sp1  <- rownames(corr_r)
corr_r_long <- gather(corr_r, sp2, corr, all_of(colnames(as.data.frame(corr_df$r))))

# Make uniq_ids
corr_r_long$uniqid <- paste(corr_r_long$sp1, corr_r_long$sp2, sep = "_")

# join both dataframes:
corr_net <- left_join(corr_r_long, corr_p_long, by = "uniqid") %>% 
  as_tibble() %>% 
  select(sp1, sp2, corr, prob)

# remove corr = 1 (correlation to same asv)
corr_net <- corr_net %>% subset(!corr == 1) 

# save nodes information
nodes <- data.frame(label = rownames(df_trans))
nodes <- left_join(nodes, nw, by = c("label" = "label")) %>% as_tibble()

############################################################################
# 4 - Keep only the positive and sup to corr_score correlation scores
############################################################################

corr_score   <- corr_threshold # minimum correlation score we use in the manuscript
corr_net_pos <- corr_net %>% subset(corr > 0) %>% subset(corr >= corr_score)

#Save the correlation results (just in case and time saving)
links_pos <- data.frame(from = corr_net_pos$sp1, 
                        to = corr_net_pos$sp2, 
                        corr = corr_net_pos$corr)

##### start ####

species_obs <- read_csv("C:/Users/roschw001/Documents/R/SDM/SDM/data_nets/species_obs.csv")


so <- species_obs[1:100,]
net <- graph_from_data_frame(so, directed = FALSE)
net_pos <- net # 

net_pos <- igraph::graph_from_data_frame(d=links_pos, directed=F) #vertices=nodes

cl <- cluster_louvain(net_pos)#, resolution = 1 # alternative used in the manuscript

# check number of communities
length(cl) # 87 communities detected

# check membership
membership(cl)

############################################################################
# 7 - Create a dataframe with membership info 
############################################################################

membercl <- data.frame(group = cl$membership, label = cl$names)

#join with species names
colnames(so)[1] <- "label"
so$label = as.character(so$label)
membercl <- membercl %>% group_by(group) %>% mutate(nb_in_group = n_distinct(label)) %>% 
  left_join(so, by = c("label" = "label"))

communities0 <- as.data.frame(cbind(membercl$group, membercl$species))
communities <- communities0[!duplicated(communities0$V1),]


###########

#write.csv(species_obs, "C:/Users/roschw001/Documents/R/SDM/SDM/data nets.csv")
# Explore the network
#V(net)
#E(net)

# Visualize the network
plot(net)

net2 <- simplify(net)

plot(net2, vertex.label = NA)

plot(net, edge.color="gray", vertex.color="lightblue", vertex.label = NA) 


#### community detection ####

ceb <- cluster_edge_betweenness(net) 
ceb2 <- cluster_walktrap(net)
ceb3 <- cluster_fast_greedy(net)
ceb4 <- cluster_leading_eigen(net)

#dendPlot(ceb, mode="hclust")

plot(ceb, net, vertex.label=NA)
plot(ceb2, net, vertex.label=NA)
plot(ceb3, net, vertex.label=NA)
plot(ceb4, net, vertex.label=NA)

ceb
ceb2
ceb3
ceb4

#edge.betweeness.community(net)

#ceb$bridges
#class(ceb)

## [1] "communities"

length(ceb)     # number of communities

## [1] 5

#membership(ceb) # community membership for each node

## s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 s11 s12 s13 s14 s15 s16 s17 

##   1   2   3   4   1   4   3   3   5   5   4   4   4   4   1   4   4

modularity(ceb) # how modular the graph partitioning is

## [1] 0.292476

#crossing(ceb, net)   # boolean vector: TRUE for edges across communities


