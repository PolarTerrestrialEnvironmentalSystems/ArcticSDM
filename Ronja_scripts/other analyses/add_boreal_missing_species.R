#### boreal change ####
#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(psych)
library(igraph)
library(ggplot2)
library(glue)
library(stringr)

#read data

getwd()
membercl1 <- read.csv("membercl1.csv")
membercl4 <- read.csv("membercl4.csv")

#filter boreal
boreal1 <- membercl1 %>% filter(str_detect(new_group, "taiga"))
boreal4 <- membercl4 %>% filter(str_detect(new_group, "taiga"))

#join, select new cells
joined_data <- boreal4 %>%
  left_join(boreal1, by = "label", suffix = c(".boreal4", ".boreal1")) %>%
  mutate(missing_in_boreal1 = is.na(new_group.boreal1)) %>%
  select(label, new_group.boreal4, new_group.boreal1, missing_in_boreal1)

summary(as.factor(joined_data$missing_in_boreal1))
#FALSE  TRUE 
#23031  4391 

new <- subset(joined_data, joined_data$missing_in_boreal1 == TRUE)

#occurence data
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"
load(glue::glue("{data}/sigma1Array.rda")) #array
sp1 <- as.data.frame(t(spArray[,,1,1]))
sp1$label = 1:38657

#species
path <- "//smb.isipd.dmawi.de/projects/bioing/user/slisovsk/ArcticSDM/Arrays_dispersalVector_Sigmoid"
files <- list.files(path, full.names = F, recursive = F)
spTable <- as.data.frame(files)
colnames(spTable) <- "species"

#old boreal forest
occ1 <- boreal1 %>%
  left_join(sp1, by = "label")
occ1b <- as.data.frame(occ1[,4:1177])
colnames(occ1b) <- spTable$species

#find out which species...
#are there
occ1c <- occ1b[, colSums(occ1b) != 0]
#are missing
occ10 <- occ1b[, colSums(occ1b) == 0] #0 species missing in the old boreal forest

occ10 <- occ1b[, colSums(occ1b) < 500]



#new boreal forest
sp4 <- as.data.frame(t(spArray[,,4,3]))
sp4$label = 1:38657

#join occs with new boreal cells
occ4 <- new %>%
  left_join(sp4, by = "label")

occ4b <- as.data.frame(occ4[,5:1178])
colnames(occ4b) <- spTable$species

#find out which species...
#are there
occ4c <- occ4b[, colSums(occ4b) != 0]
#are missing
occ40 <- occ4b[, colSums(occ4b) == 0] #71 species missing in new boreal forest


s <- as.data.frame(colnames(occ10))
colnames(s) <- "species"
misssp <- as.data.frame(colnames(occ40))
colnames(misssp) <- "species"
m <- misssp %>%
  left_join(s, by = "species")


#traits
setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM")
path <- "C:/Users/roschw001/Documents/migration values"
traits <- read.csv("C:/Users/roschw001/Documents/migration values/migration_meters_all.csv", dec=",")
traits$height <- gsub(",", ".", traits$height)
traits$height <- as.numeric(traits$height)
h0 <- read.csv("hypotheses_results/H0_new/results_sigma1new.csv", sep=";")
h0t <- h0 %>%
  left_join(traits, by = "species")
h0t$distance99[134] <- 1500
h0t$gform[134] <- "Shrub"
h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"



spt <- misssp %>%
  left_join(h0t, by = "species") 

summary(spt)

summary(as.factor(spt$gform))
summary(as.factor(spt$mode))
summary(as.factor(spt$class))
summary(spt$height)
summary(as.factor(spt$family))
summary(as.numeric(spt$distance99))

subset(spt, spt$class == "4")
