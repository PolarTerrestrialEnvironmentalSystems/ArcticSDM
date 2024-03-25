#mars test
#https://uc-r.github.io/mars model

#### input ####
#for model
envtrain <- readRDS("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/envtrain.Rds")

#for evaluation
testpres <- readRDS("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/testpres.Rds")
testbackg <- readRDS("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/testbackg.Rds")

#### model ####
library(earth)

mars1 <- earth(
  pa ~ .,  
  data = envtrain)   

#### evalution ####

print(mars1)
evaluate(testpres, testbackg, mars1)
