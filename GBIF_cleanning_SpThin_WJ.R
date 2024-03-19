# Cleaning GBIF data
# Original codes: https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_GBIF_data_with_CoordinateCleaner.html

#----------------------------------------------------------------------------
library(devtools)
#install_github("ropensci/CoordinateCleaner")
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sf)
library(rnaturalearthdata)

# import the occurrence data downloaded from GBIF
gbif <- read.csv("E:/AWI/glacial_refugia_project/Framework/GBIF/Betula_0081065-231120084113126.csv",
                 sep="\t",fill=TRUE,header=TRUE, quote="",encoding="UTF-8")

# check the species list
gbif %>% group_by(species) %>% summarise(count = n_distinct(gbifID))

# select columns of interest
dat <- gbif %>%
  dplyr::select(species, decimalLongitude, 
                decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters,
                year, basisOfRecord, institutionCode)

# remove records without coordinates
dat <- dat %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

# convert country code from ISO2c to ISO3c
dat$countryCode <-  countrycode(dat$countryCode, 
                                origin =  'iso2c',
                                destination = 'iso3c')

# remove flag errors
# flag errors including: sea coordinates, zero coordinates,
# coordinate - country mismatches, coordinates assigned to country and province centroids,
# coordinates within city areas, outlier coordinates and coordinates assigned to
# biodiversity institutions
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids",
                                     "equal", "zeros", "countries")) # most test are on by default

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude problematic records
dat_cl <- dat[flags$.summary,]

# the flagged records
dat_fl <- dat[!flags$.summary,]

# remove records with a precision below 1 km
# first get an overview
count <- dat_cl %>% group_by(coordinateUncertaintyInMeters) %>% summarise(count = n_distinct(gbifID))

dat_cl %>% 
  mutate(Uncertainty = coordinateUncertaintyInMeters / 1000) %>% 
  ggplot(aes(x = Uncertainty)) + 
  geom_histogram() +
  xlab("Coordinate uncertainty in meters") +
  theme_bw()

dat_cl <- dat_cl %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 1 | is.na(coordinateUncertaintyInMeters))

# so what about the data without the information 'coordinateUncertaintyInMeters'?
# function to check if the coordinates after the decimal point are less than two digits (must be >0)
check_pattern <- function(coord) {
  grepl("\\.\\d{2,}$", as.character(coord))
}

# Filter rows based on the pattern
dat_cl_filter <- subset(dat_cl, check_pattern(dat_cl$decimalLongitude) & check_pattern(dat_cl$decimalLatitude))

dat_cl_removed <- subset(dat_cl, !(check_pattern(dat_cl$decimalLongitude) & check_pattern(dat_cl$decimalLatitude)))


# remove unsuitable data sources, especially fossils 
table(dat_cl_filter$basisOfRecord)

dat_cl_filter <- filter(dat_cl_filter, basisOfRecord != "FOSSIL_SPECIMEN")

# individual count > 0
table(dat_cl_filter$individualCount)

dat_cl_filter <- dat_cl_filter %>%
  filter(individualCount > 0 | is.na(individualCount)) # high counts are not a problem

# age of records
table(dat_cl_filter$year)

dat_cl_filter <- dat_cl_filter %>% filter(year >= 1970)


# check the family and taxon rank
table(dat_cl_filter$family) # only the family you want?
table(dat_cl_filter$taxonRank)


# check the cleaned species list
taxa_list <- dat_cl_filter %>% group_by(species) %>% summarise(count = n_distinct(gbifID))

# remove species occurrences less than 10
dat_cl_filter_coo10 <- subset(dat_cl_filter, !(species %in% c("Betula alaskana", "Betula aurata",
                                                              "Betula avatshensis", "Betula dugleana",
                                                              "Betula intermedia", "Betula paramushirensis",
                                                              "Betula pendula × pubescens","Betula pumila",
                                                              "Betula utahensis")))

# check again
dat_cl_filter_coo10 %>% group_by(species) %>% summarise(count = n_distinct(gbifID))


# select the columns we need
dat_cl_filter_coordinates <- dat_cl_filter_coo10 %>%
  dplyr::select(gbifID, decimalLongitude, decimalLatitude)

# rename the columns
colnames(dat_cl_filter_coordinates)[2] <- 'Longitude'
colnames(dat_cl_filter_coordinates)[3] <- 'Latitude'

# remove duplicates
occurence <- unique(dat_cl_filter_coordinates[, c("Longitude", "Latitude")])

# plot the cleaned data on map
wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() + coord_fixed() + wm +
  geom_point(data = occurence,
             aes(x = Longitude, y = Latitude),
             colour = "darkred",
             size = 0.5) + theme_bw()

# save the cleaned data
write.table(occurence, file="E:/AWI/glacial_refugia_project/Framework/GBIF/gbif_betula_cleaned_OnlyCoordinates_29.01.2024.csv",
            append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# Spatial thinning
# Original codes: https://cran.r-project.org/web/packages/spThin/vignettes/spThin_vignette.html

#----------------------------------------------------------------------------
library(spThin)

# create a new column named by the genus/family name
occurence$SPEC <- 'Betula'

# spatial thinning to 1 km
thinned_dataset_full <-
  thin( loc.data = occurence, 
        lat.col = "Latitude", long.col = "Longitude", 
        spec.col = "SPEC", 
        thin.par = 1, #the distance (in km) that you want records to be separated by
        reps = 100, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 1, 
        out.dir = "E:/AWI/glacial_refugia_project/Framework/GBIF/", out.base = "betula_thinned", 
        write.log.file = TRUE,
        log.file = "thinned_full_log_file.txt" )

plotThin(thinned_dataset_full)

# plot the thinned data on map
wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() + coord_fixed() + wm +
  geom_point(data = thin,
             aes(x = LONG, y = LAT),
             colour = "darkred",
             size = 0.00001) +
  theme_bw()

# save the thinned data
write.table(thinned_dataset_full, file="E:/AWI/glacial_refugia_project/Framework/GBIF/gbif_betula_cleaned_thinned_29.01.2024.csv",
            append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)