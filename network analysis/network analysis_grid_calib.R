#network analysis

#### packages ####
library(stars)
library(tidyr)
library(terra)
library(rgbif)
library(sf)
sf_use_s2(FALSE)
library(tidyverse)

#### input ####
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM"
path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/Results"

data <- "D:"
path <- "D:/Results"

data <- "C:/Users/roschw001/Documents/R/SDM/SDM/festplatte_sub"
#Achillea millefolium_MaxEnt_calibration.tif"

#region
ecoreg <- st_read(glue::glue("{data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
  filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0)

reg <- subset(ecoreg, ecoreg$ECO_ID_U == "17077") #choose alaska



#### species names ####
#list.dirs
flsl      <- list.files(data, 
                       recursive = FALSE, full.names = T)

speciesnames = as.data.frame(matrix(nrow = 1, ncol=1))

# Split the path by the "/" character
for (i in 1:length(flsl)){

# Extract the last part of the path, which contains the species names
path_parts <- unlist(strsplit(flsl[i], "/|_"))

# Extract the last part of the path, which contains the species name
speciesnames[i,1] <- path_parts[10] #[3]

}


colnames(speciesnames)[1] <- "species"


#for(species in speciesnames$species){

#### read & grid & cor data ####
grids <- list()

for(sp in 1:3){

species = speciesnames[sp,1]

dat <- read_stars(glue::glue("{data}/Results/{species}_MaxEnt_calibration.tif"))  

grids[[sp]] <- dat %>% st_make_grid(., cellsize = 100000) %>%
            st_extract(dat, ., cellsize = 100000) %>%
            st_as_sf(.) %>% st_drop_geometry()
    
}

tab <- do.call(cbind, grids)
tab <- replace(tab, is.na(tab), 0)

cormatrix <- cor(tab)
cormatrix

#### create net ####
library(igraph)
net_all <- graph_from_data_frame(cormatrix, directed = FALSE)
plot(net_all)
#### clustering ####
cl <- cluster_louvain(net_all) 

#### output ####.
# check number of communities
length(cl) #95 communities detected
cl
# check membership
membership(cl)

membercl <- data.frame(group = cl$membership, label = cl$names)
