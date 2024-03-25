#glm test

#### input ####
library(dismo)

#for model
envtrain <- readRDS("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/envtrain.Rds")

#for evaluation
testpres <- readRDS("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/testpres.Rds")
testbackg <- readRDS("//smb.isipd.dmawi.de/projects/bioing/data/ArcticSDM/test/testbackg.Rds")


#### predict ####

# #with raster data
# predictors <- stack(list.files(file.path(system.file(package="dismo"), 'ex'), pattern='grd$', full.names=TRUE ))
# p <- predict(predictors, m1)
# plot(p)

#### model with evaluation ####

gm1 <- glm(pa ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17,
           family = binomial(link = "logit"), data=envtrain)

#### evaluation ####
summary(gm1) 
coef(gm1)
evaluate(testpres, testbackg, gm1)



