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
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/"

#load(glue::glue("{data}/spTable.rda"))
#load(glue::glue("{data}/spArray_5km.rda"))

load(glue::glue("{data}/sdm_occurance_array_array_25km.rda")) #only array
load(glue::glue("{data}/sdm_occurance_array_metadata_25km.rda")) #only meta

spTable <- as.data.frame(spArrayMeta[[1]])
sppres <- as.data.frame(spArray[1,,,])

#### richnesstab ####

richnesstab <-
  as.data.frame(matrix(nrow = ncol(sppres), ncol=12))
colnames(richnesstab) <- c("y1s1", "y2s1", "y3s1", "y4s1", 
                          "y1s2","y2s2","y3s2", "y4s2",
                          "y1s3","y2s3","y3s3","y4s3")

#s1
for (y in 1:4){
sppres <- as.data.frame(spArray[,,y,1])
  richnesstab[,y] <- as.data.frame(apply(sppres, 2, sum))
}


#s2
for (y in 1:4){
  sppres <- as.data.frame(spArray[,,y,2])
  richnesstab[,y+4] <- as.data.frame(apply(sppres, 2, sum))
}

#s3
for (y in 1:4){
sppres <- as.data.frame(spArray[,,y,3])
richnesstab[,y+8] <- as.data.frame(apply(sppres, 2, sum))
  }


#### summary ####

library(psych)
library(tibble)

# Example dataframe

# Calculate summary statistics and save as a new tibble
y1s1 <- describe(richnesstab$y1s1) %>%
  as_tibble(rownames = "Variable")

# Calculate summary statistics for each column
summary_stats <- lapply(richnesstab, function(x) {
  describe(x) %>%
    as_tibble(rownames = "Variable")
})

# Combine summary statistics into a single tibble
summary_tbl <- bind_rows(summary_stats, .id = "Column")

summary_tbl$year <- c(1,2,3,4,1,2,3,4,1,2,3,4)
summary_tbl$yearreal <- c(1996, 2026,2056,2086,1996, 2026,2056,2086,1996, 2026,2056,2086)
summary_tbl$yearend <- c(2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100, 2010, 2040, 2070, 2100)
summary_tbl$scen <- c(1,1,1,1,2,2,2,2,3,3,3,3)
summary_tbl$scenario <- c("ssp126", "ssp126","ssp126","ssp126","ssp370", "ssp370","ssp370","ssp370","ssp585", "ssp585", "ssp585","ssp585")

summary_tbl$scen <- as.factor(summary_tbl$scen)
summary_tbl$scenario <- as.factor(summary_tbl$scenario)
getwd()

write.csv2(summary_tbl, "H2_summary_tbl.csv")
##### mean ####

mod_y <- lm(summary_tbl$mean ~ summary_tbl$year)
summary(mod_y)
plot(summary_tbl$year, summary_tbl$mean)
abline(mod_y, col="red")

library(ggplot2)
library(dplyr)


# Create the plot
ggplot(summary_tbl, aes(yearend, mean)) +
  geom_point(aes(color = scenario)) +
  geom_smooth(aes(color = scenario), method = "lm", se = FALSE) +
  scale_color_manual(values = c("ssp126" = "blue", "ssp370" = "orange", "ssp585" = "red")) +
  theme_minimal() 

plot(summary_tbl$scen, summary_tbl$mean)
mod_s <- lm(summary_tbl$mean ~ summary_tbl$scen)
summary(mod_s)

##### max ####

mod_y <- lm(summary_tbl$max ~ summary_tbl$year)
summary(mod_y)
plot(summary_tbl$year, summary_tbl$max)
abline(mod_y, col="red")

plot(summary_tbl$scen, summary_tbl$max)
mod_s <- lm(summary_tbl$max ~ summary_tbl$scen)
summary(mod_s)

ggplot(summary_tbl, aes(yearend, max)) +
  geom_point(aes(color = scenario)) +
  geom_smooth(aes(color = scenario), method = "lm", se = FALSE) +
  scale_color_manual(values = c("ssp126" = "blue", "ssp370" = "orange", "ssp585" = "red")) +
  theme_minimal() 
warnings()

