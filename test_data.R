#testdata for models

#### input ####
##### env #####
path <- file.path(system.file(package="dismo"), 'ex')
library(dismo)
files <- list.files(path, pattern='grd$', full.names=TRUE)
predictors <- stack(files)

##### occ #####
library(maptools)
data(wrld_simpl)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file,  header=TRUE,  sep=',')
# we do not need the first column
bradypus  <- bradypus[,-1]

#### data train & test ####
#extent
ext <- extent(-90, -32, -33, 23)

#test & training split
set.seed(0)
group <- kfold(bradypus, 5)
pres_train <- bradypus[group != 1, ]
pres_test <- bradypus[group == 1, ]

#background points
set.seed(10)
pred_nf <- dropLayer(predictors, 'biome')
backg <- randomPoints(pred_nf, n=1000, ext=ext, extf = 1.25)
colnames(backg) = c('lon', 'lat')
group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

#create data set
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(predictors, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain[,'biome'] = factor(envtrain[,'biome'], levels=1:14)

testpres <- data.frame( extract(predictors, pres_test) )
testbackg <- data.frame( extract(predictors, backg_test) )
testpres[ ,'biome'] = factor(testpres[ ,'biome'], levels=1:14)
testbackg[ ,'biome'] = factor(testbackg[ ,'biome'], levels=1:14)

#### export ####
saveRDS(envtrain, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/envtrain.Rds")
saveRDS(testpres, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/testpres.Rds")
saveRDS(testbackg, "//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/testbackg.Rds")