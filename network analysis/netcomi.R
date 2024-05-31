#netcomi


#### installation ####
#https://github.com/stefpeschel/NetCoMi?tab=readme-ov-file

rm(list=ls()) #clear
# 
# # Required packages
# install.packages("devtools")
# install.packages("BiocManager")
# 
# 
# # Install NetCoMi
# devtools::install_github("stefpeschel/NetCoMi", 
#                          dependencies = c("Depends", "Imports", "LinkingTo"),
#                          repos = c("https://cloud.r-project.org/",
#                                    BiocManager::repositories()))
# 
# 
# unloadNamespace("miniUI")
# unloadNamespace("shiny")
# unloadNamespace("htmltools")
# 
# devtools::install_github("zdk123/SpiecEasi")
# devtools::install_github("GraceYoon/SPRING")
# 
# installNetCoMiPacks()

library(NetCoMi)

library(tidyr)
library(tidyverse)

#### data ####
#test data
library(NetCoMi)
data("amgut1.filt")
data("amgut2.filt.phy")

amgut1.filt

library(stars)
load("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM_data/GBIF/GBIF_dataComp/occTab_all.rda")
head(occTab_all)

species <- occTab_all$species
coords <- occTab_all$geometry


coordsxy <- as.data.frame(st_coordinates(occTab_all, dims = c("x", "y")))

summary(coordsxy)
coordsxy$x1 = round(coordsxy$X, -4)
coordsxy$y1 = round(coordsxy$Y, -4)
nettab = cbind(coordsxy[,3:4], species)


# Unite 'col1' and 'col2' into a single column 'united_col'
nettab <- unite(nettab, col = site, x1, y1, sep = "_")
nettab$p = 1

#nettab$site = as.numeric(rownames(nettab))
#head(nettab2)
#nettab1 = unlist(nettab)


#a <- nettab$species
# Transform the DataFrame to wide format

#nettab3 <- nettab %>% group_by(across(species)) 
  
#nettabl <- nettab %>%
 # pivot_wider(names_from = species, values_from = p)

nettab[1:10,1]

#nettab$s <- as.factor(nettab$site)


nrow(nettab)
657314
# Assuming you have a long-format DataFrame 'df'
df <- data.frame(id = nettab[1:10,1],
                 year = c("Abies balsamea", "Abies balsamea", "Abies balsamea","Abies balsamea",
                          "Abies balsamea","Abies balsamea","Abies balsamea","Abies balsamea", 
                          "Abietinella abietina", "Abutilon theophrasti"),
                 value = c(1,1,1,1,1,1,1,1,1,1))


#### subset data ####
library(dplyr)
df = nettab 
df$site
df3 <- df %>% filter(site == "4830000_140000")
df4 <- df %>% filter(site %in% c("4830000_140000", "4900000_160000", "4830000_50000",
                                 "4830000_110000", "4890000_140000", "4890000_110000", "4870000_90000", 
                                 "4850000_70000",  "4910000_160000", "4830000_120000", "4850000_60000",  
                                 "4920000_160000", "4830000_110000", "4910000_1e+05"))
unique(df$species)


dfs <- df %>% filter(species %in% c("Picea abies",                   "Picea glauca",                  "Picea mariana",                
                                 "Picea pungens",                 "Picea rubens",                  "Pilea pumila",                 
                                  "Pilosella aurantiaca",          "Pilosella caespitosa",          "Pilosella officinarum",        
                                  "Pilosella piloselloides",       "Pinus banksiana",               "Pinus resinosa",               
                                  "Pinus strobus",                 "Pinus sylvestris",              "Piptatheropsis pungens",       
                                  "Plagiochila porelloides",       "Plagiomnium ciliare",           "Plagiomnium cuspidatum",       
                                 "Plagiothecium cavifolium",      "Plagiothecium denticulatum",    "Plagiothecium laetum"))

a = as.data.frame(table(df))


a1 = sort(a, a$Freq)
df <- data.frame(site = nettab[1:10000,1],
                 species = nettab[1:10000,2],
                 p = nettab[1:10000,3])

#nettab4 <- nettab %>% filter()
  
# Transform the DataFrame to wide format
netwide <- df %>%
  pivot_wider(names_from = species, values_from = p,
              values_fn = list(p = sum)) 



netwide0 <- replace(netwide, is.na(netwide), 0)

#summary(netwide0)

netwide1 = netwide0[1:1000, 1:100]
netwide2 = netwide0[1:1000, 2:100]
netwide1 <- na.omit(netwide0)
#construct net
net_spring <- netConstruct1(netwide1,
                           filtTax = "highestFreq", # "none",
                           filtTaxPar = list(highestFreq = 10),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1),
                           measure = "spieceasi", #spring,
                           measurePar = NULL,
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = NULL,
                           jointPrepro = FALSE,)


wnet_spring2 <- netConstruct(netwide1,
                            data2 = NULL,
                            #dataType = "counts",
                            #group = NULL,
                            #matchDesign = NULL,
                           # taxRank = NULL,
                            
                            # Association/dissimilarity measure:
                            measure = "spieceasi",
                           # measurePar = NULL,
                            
                            # Preprocessing:
                            jointPrepro = TRUE,
                            filtTax = "none",
                            filtTaxPar = NULL,
                            filtSamp = "none",
                            filtSampPar = NULL,
                            zeroMethod = "none",
                            zeroPar = NULL,
                            normMethod = "none",
                            normPar = NULL,
                            
                            # Sparsification:
                            sparsMethod = "t-test",
                            thresh = 0.3,
                            alpha = 0.05,
                            adjust = "adaptBH",
                            trueNullMethod = "convest",
                            lfdrThresh = 0.2,
                            nboot = 1000L,
                            assoBoot = NULL,
                            cores = 1L,
                            logFile = "log.txt",
                            softThreshType = "signed",
                            softThreshPower = NULL,
                            softThreshCut = 0.8,
                            kNeighbor = 3L,
                            knnMutual = FALSE,
                            
                            # Transformation:
                            dissFunc = "signed",
                            dissFuncPar = NULL,
                            simFunc = NULL,
                            simFuncPar = NULL,
                            scaleDiss = TRUE,
                            weighted = TRUE,
                            
                            # Further arguments:
                            sampleSize = NULL,
                            verbose = 2,
                            seed = NULL)
#analyse cor-Matrix
props_spring <- netAnalyze(net_spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

plotHeat(mat = props_spring$graphletLCC$gcm1,
         pmat = props_spring$graphletLCC$pAdjust1,
         type = "mixed",
         title = "GCM", 
         colorLim = c(-1, 1),
         mar = c(2, 0, 2, 0))

# Add rectangles highlighting the four types of orbits
graphics::rect(xleft   = c( 0.5,  1.5, 4.5,  7.5),
               ybottom = c(11.5,  7.5, 4.5,  0.5),
               xright  = c( 1.5,  4.5, 7.5, 11.5),
               ytop    = c(10.5, 10.5, 7.5,  4.5),
               lwd = 2, xpd = NA)

text(6, -0.2, xpd = NA, 
     "Significance codes:  ***: 0.001;  **: 0.01;  *: 0.05")


#?summary.microNetProps
summary(props_spring, numbNodes = 5L)

# help page
?plot.microNetProps

p <- plot(props_spring, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network on OTU level with SPRING associations", 
          showTitle = TRUE,
          cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)

p$q1$Arguments$cut

#### Export to Gephi ####

# For Gephi, we have to generate an edge list with IDs.
# The corresponding labels (and also further node features) are stored as node list.

# Create edge object from the edge list exported by netConstruct()
edges <- dplyr::select(net_spring$edgelist1, v1, v2)

# Add Source and Target variables (as IDs)
edges$Source <- as.numeric(factor(edges$v1))
edges$Target <- as.numeric(factor(edges$v2))
edges$Type <- "Undirected"
edges$Weight <- net_spring$edgelist1$adja

nodes <- unique(edges[,c('v1','Source')])
colnames(nodes) <- c("Label", "Id")

# Add category with clusters (can be used as node colors in Gephi)
nodes$Category <- props_spring$clustering$clust1[nodes$Label]

edges <- dplyr::select(edges, Source, Target, Type, Weight)

setwd("C:/Users/roschw001/Documents/R/SDM/SDM/netcomi")
write.csv(nodes, file = "nodes.csv", row.names = FALSE)
write.csv(edges, file = "edges.csv", row.names = FALSE)
