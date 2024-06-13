
#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)

#input
setwd("//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays")
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/sdm_arrays"

load(glue::glue("{data}/spTable.rda"))
load(glue::glue("{data}/binaryArray_25km.rda"))

#binarArray[species, cells, year, scenario]
sp1_pres <- binarArray[1,,1,1]
sp1_fut  <- binarArray[1,,4,3]

sp1a <- binarArray[1,,,]
resultlist <- list()

for(sp in 1:15){

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

resultlist[[sp]] <- results
}

resultlist[[15]]

sp1
summary(sp1_pres)
summary(sp1_fut)
sum(sp1_pres)
sum(sp1_fut)

apply(sp1, 1:2, sum)
