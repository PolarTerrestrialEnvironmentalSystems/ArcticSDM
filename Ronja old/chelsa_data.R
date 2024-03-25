#data ulrike
library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)
wd <- "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/data/chelsa_data_links"


fls      <- tibble(fls = list.files(wd,
                                    pattern = ".txt$",  
                                    recursive = FALSE))

for(t in fls){
  model_list <- read_stars(glue::glue("{wd}/{t}"))}

# Create a vector of URLs

urls = read.table(glue::glue("{wd}/envidatS3paths_2011_2040_126.txt"))
head(url)
# Create a directory to store the downloaded files
dir.create(glue::glue("{wd}/chelsatrace21k"))

# Loop through the URLs and download the files
for (url in urls[1]) {
  filename <- basename(url)
  download.file(url, destfile = paste0("chelsatrace21k/", filename))
}


