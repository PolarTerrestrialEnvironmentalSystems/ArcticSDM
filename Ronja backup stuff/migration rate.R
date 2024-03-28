#migration rate

#### read data ####
rates_zani <- read.csv("C:/Users/roschw001/Documents/migration values/migrationrates_zani2023.csv", sep=";", dec=".")
rates_meier <- read.csv("C:/Users/roschw001/Documents/migration values/migrationrates_meier2011.csv", sep=";", dec=".")

traits <- read.csv("C:/Users/roschw001/Documents/migration values/sibala_traits_google_height.csv - sibala_traits_summary.csv", sep=",", dec=",")
dispersal <- read.csv("C:/Users/roschw001/Documents/migration values/sibala_traits_google_height.csv - dispersal_Lososova.csv", sep=",", dec=",")

species_mod <- read.csv("C:/Users/roschw001/Documents/migration values/SpeciesList.csv", sep=",", dec=".")

species_mod50 = species_mod %>% filter(nrOcc > 1000)

speciestraits <- species_mod50 %>% left_join(traits, by = join_by("species" == "species"))
speciestraitsna = na.omit(speciestraits) #we have data

#50  3221 1262
#100 2381 1111
#200 1675  919
#300 1304  775
#500  899  613
#1000 494  400

summary(species_mod$nrOcc)
#### join ####
library(tidyverse)
rates0 <- rates_meier %>% full_join(rates_zani, by = join_by("species" == "species"), 
                                   relationship= "one-to-many") 

#### try to extract unique species--> failed ####
rates0 = rates0 %>% arrange(species)
species = as.data.frame(unique(rates0$species)) #26
unique(species)
colnames(species)[1] = "speciess"
species
species = species %>% arrange(speciess)
species2 = names(table(species$species))
species$speciess = as.character(species$speciess)
unique_values <- species %>% distinct(speciess) %>% pull(speciess)


#### join seperately ####
#rates1 <- species %>% left_join(rates_meier, by = join_by("species" == "species"))
#rates2 <- rates1 %>% left_join(rates_zani, by = join_by("species" == "species"))

ratestraits <- rates_meier %>% left_join(traits, by = join_by("species" == "species"))
ratestraitsz <- rates_zani %>% left_join(traits, by = join_by("species" == "species"))

#### select relevant columns and rows ####
rtz = ratestraitsz %>% select(c(species, mig_km_y_free, Height..m....all.data.with.SearchBasedData...BET,
                                Seed.mass..mg....all.data, Growth.form))
  
rtz3 = rtz %>% filter(species %in% c("Picea abies","Ulmus glabra"))
colnames(rtz3) = c("species", "mean", "height", "seed.mass", "growth.form")  
                              
rt = ratestraits[3:7,]
rt3 = rt   %>% select(c(species, mean, Height..m....all.data.with.SearchBasedData...BET,
                        Seed.mass..mg....all.data, Growth.form))
colnames(rt3) = c("species", "mean", "height", "seed.mass", "growth.form")  

rt3 = rt3 %>% arrange(species)
rt4 = rt3[2:5,]

rtboth =  bind_rows(rt4, rtz3)
rtboth = rtboth %>% arrange(species)
rtboth$mean[3] <- (10.5+3.2)/2 #only once picea abies, mig average

rt_data = rtboth[-4,]
rt_data = rt_data %>% arrange(species)

#### join with dispersal ####
rt_datad <- rt_data %>% left_join(dispersal, by = join_by("species" == "species"))

rt_vars = rtboth[-4,2:4] #only vars

rt_datad2 = cbind(rt_datad[,1:5], rt_datad[,9:13])

#### cor analysis ####
library(ellipse)
            
cc <- cor(rt_vars, method="spearman")
plotcorr(cc, cex.lab=0.5)

cor.test(rt_vars$mean, rt_vars$height)
cor.test(rt_vars$mean, rt_vars$seed.mass)




