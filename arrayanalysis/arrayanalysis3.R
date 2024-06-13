#1 with table 12 rows
#2 with table 4 rows, per scenario
#3 with species list and tab, 12 rows
#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)

#input
setwd("//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays")
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays"

load(glue::glue("{data}/spTable.rda"))
load(glue::glue("{data}/binaryArray_25km.rda"))

#binarArray[species, cells, year, scenario]
sp1_pres <- binarArray[1,,1,1]
sp1_fut  <- binarArray[1,,4,3]

sp1a <- as.data.frame(binarArray[1,,,])





resultlist <- list()
percenttab <- as.data.frame(matrix(nrow = nrow(spTable), ncol=14))
colnames(percenttab) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                       "y1s2","y2s2","y3s2", "y4s2",
                       "y1s3","y2s3","y3s3","y4s3", "trend", "species")

areatab <- as.data.frame(matrix(nrow = nrow(spTable), ncol=13))
colnames(areatab) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                          "y1s2","y2s2","y3s2", "y4s2",
                          "y1s3","y2s3","y3s3","y4s3", "species")

#### base: areatab & percenttab ####
for(sp in 1:nrow(spTable)){


sp1a <- as.data.frame(binarArray[sp,,,])

results <- as.data.frame(matrix(ncol=12, nrow=3))

colnames(results) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                       "y1s2","y2s2","y3s2", "y4s2",
                       "y1s3","y2s3","y3s3","y4s3")

rownames(results) <- c("sum", "delta", "percent")


#sum
for (i in 1:ncol(results)) {
  # Calculate the sum of values in the first 38657 rows for the current column
  sum_val <- sum(sp1a[1:38657, i])
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[1, i] <- sum_val
}



#delta
for (i in 1:3) {
  # Calculate the sum of values in the first 38657 rows for the current column
  delta_val <- results[1, i+1] - results[1, 1]
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[2, i+1] <- delta_val
}

for (i in 1:3) {
  # Calculate the sum of values in the first 38657 rows for the current column
  delta_val <- results[1, i+5] - results[1, 1]
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[2, i+5] <- delta_val
}

for (i in 1:3) {
  # Calculate the sum of values in the first 38657 rows for the current column
  delta_val <- results[1, i+9] - results[1, 1]
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[2, i+9] <- delta_val
}

#percent
for (i in 1:ncol(results)) {
  # Calculate the sum of values in the first 38657 rows for the current column
  per_val <- results[2,i]/results[1,1]*100
  
  # Assign the sum to the 38658th row (the last row) of the current column
  results[3, i] <- per_val
}

#resultlist
resultlist[[sp]] <- results

#### percenttab H0 ####
percenttab[sp,] <- results[3,]
percenttab[sp,14] <- spTable[sp, 1]

percenttab[sp,13] <- ifelse(mean(results[3,4], results[3,8], results[3,12]) > 0, "increase",
                        ifelse(mean(results[3,4], results[3,8], results[3,12]) < 0, "decrease",
                               ifelse(mean(results[3,4], results[3,8], results[3,12]) == 0, "stable",
                                      NA)))
percenttab[sp, 15] <- ifelse(results[3,4]> 0, "increase",
                             ifelse(results[3,4]< 0, "decrease",
                                    ifelse(results[3,4] == 0, "stable",
                                           NA)))
percenttab[sp, 16] <- ifelse(results[3,8]> 0, "increase",
                             ifelse(results[3,8]< 0, "decrease",
                                    ifelse(results[3,8] == 0, "stable",
                                           NA)))
percenttab[sp, 17] <- ifelse(results[3,12]> 0, "increase",
                             ifelse(results[3,12]< 0, "decrease",
                                    ifelse(results[3,12] == 0, "stable",
                                           NA)))
colnames(percenttab) [15:17] <- c("trend_s1", "trend_s2", "trend_s3")
#areatab H1
areatab[sp,] <- results[1,]
areatab[sp,13] <- spTable[sp, 1]

}


summary(as.factor(percenttab$trend))

#### gform H1 based on areatab ####
traits <- read.csv("C:/Users/roschw001/Documents/migration values/modspeciesdist.csv", sep=",", dec=",")
traits_gform <- as.data.frame(cbind(traits$species, traits$gform))
colnames(traits_gform) <- c("species", "gform")
areatab <- areatab %>% left_join(traits_gform, by="species") 

summary(as.factor(gform_result$gform))
gform_result <- as.data.frame(matrix(nrow=7, ncol=14))  
colnames(gform_result) = colnames(areatab)
gform_result$gform <-  unique(traits_gform$gform)

##### gformarea #####
library(dplyr)

#sum area for all yxsx
r <- seq(1,13,1)
areatab <- rbind(r, areatab)

areatab[1,13] <- "sum"
for (i in 1:12){
  # Calculate the sum of valuesin the first 38657 rows for the current column
  sum_val <- sum(areatab[2:17, i])
  areatab[1, i] <- sum_val
}

gformarea <- areatab %>%
  group_by(gform) %>%
  summarise(across(1:12, ~ sum(.x, na.rm = TRUE)))
areatab %>% filter(gform)


for (i in 2:13){
gformarea[6,i] <- sum(gformarea[1:5, i])
}
gformarea[6,1] = 'all'

##### gformpercent #####
gformpercent = gformarea

for (n in 1:nrow(gformpercent)) {
  for (i in 2:13) {
    
    
    gformpercent[n,i] <- gformarea[n, i]/gformarea[6, i]*100
  }}

##### trend #####
#trend s1
for (i in 1:nrow(gformpercent)){
  gformpercent[i,14] <- ifelse(gformpercent[i,2] < gformpercent[i,5], "increase2100",
                         ifelse(gformpercent[i,2] > gformpercent[i,5], "decrease2100",
                         ifelse(gformpercent[i,2] == gformpercent[i,5], "stable", NA)))
          
}

#trend s2
for (i in 1:nrow(gformpercent)){
  gformpercent[i,15] <- ifelse(gformpercent[i,6] < gformpercent[i,9], "increase2100",
                               ifelse(gformpercent[i,6] > gformpercent[i,9], "decrease2100",
                                      ifelse(gformpercent[i,6] == gformpercent[i,9], "stable", NA)))
  
}

#trend s3
for (i in 1:nrow(gformpercent)){
  gformpercent[i,16] <- ifelse(gformpercent[i,10] < gformpercent[i,13], "increase2100",
                               ifelse(gformpercent[i,10] > gformpercent[i,13], "decrease2100",
                                      ifelse(gformpercent[i,10] == gformpercent[i,13], "stable", NA)))
  
}

colnames(gformpercent) [14:16] <- c("trend_s1", "trend_s2", "trend_s3")

summary(as.factor(gformpercent$trend_s1))
summary(as.factor(gformpercent$trend_s2))
summary(as.factor(gformpercent$trend_s3))
