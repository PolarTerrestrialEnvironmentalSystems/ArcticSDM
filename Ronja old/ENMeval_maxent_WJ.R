WD <- "E:/AWI/glacial_refugia_project/Framework/MaxEnt/test/" 

taxaname <- 'Alnus'

library(dismo)
library(maptools)
library(ENMeval)

#library(remotes)
#install_version("maptools", "0.9-5")


#=== Data preparation ====
# import environmental data
path <- file.path(system.file(package="dismo"), 'ex')
files <- list.files(path, pattern='grd$', full.names=TRUE )

predictors <- stack(files)
names(predictors)
plot(predictors)

# import occurance data
data(wrld_simpl)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file,  header=TRUE,  sep=',')
bradypus  <- bradypus[,-1]


#=== ENMevaluation ====
bg <- dismo::randomPoints(predictors[[1]], n = 10000) %>% as.data.frame()

colnames(bg) <- colnames(bradypus)

enmeval_results <- ENMeval::ENMevaluate(occs = bradypus, envs = predictors, bg = bg, 
                                        tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
                                        partitions = "randomkfold", partition.settings = list(kfolds = 2), 
                                        algorithm = "maxnet")

enmeval <- enmeval_results@results[which.min(enmeval_results@results$delta.AICc), c(1, 2)]

fc <- as.character(enmeval$fc)
rm <- as.numeric(enmeval$rm)



#========== MaxEnt ==========
set.seed(0)
group <- kfold(bradypus, 5)
pres_train <- bradypus[group != 1, ]
pres_test <- bradypus[group == 1, ]

set.seed(10)
pred_nf <- dropLayer(predictors, 'biome')
backg <- randomPoints(pred_nf, n=1000, extf = 1.25)

# maxent ready?
maxent()

# set up the argument list
args_list <- c('responsecurves=TRUE',
               'jackknife=TRUE',
               'pictures=TRUE',
               'autofeature=FALSE',
               'linear=TRUE',
               'quadratic=TRUE',
               'product=TRUE',
               'threshold=TRUE',
               'hinge=TRUE',
               paste0('betamultiplier=', rm))

# adjust arguments based on fc
if (fc == "L") {
  args_list[5:9] <- c('linear=TRUE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'hinge=FALSE')
} else if (fc == "LQ") {
  args_list[5:9] <- c('linear=TRUE', 'quadratic=TRUE', 'product=FALSE', 'threshold=FALSE', 'hinge=FALSE')
} else if (fc == "H") {
  args_list[5:9] <- c('linear=FALSE', 'quadratic=FALSE', 'product=FALSE', 'threshold=FALSE', 'hinge=TRUE')
} else if (fc == "LQH") {
  args_list[5:9] <- c('linear=TRUE', 'quadratic=TRUE', 'product=FALSE', 'threshold=FALSE', 'hinge=TRUE')
} else if (fc == "LQHP") {
  args_list[5:9] <- c('linear=TRUE', 'quadratic=TRUE', 'product=TRUE', 'threshold=FALSE', 'hinge=TRUE')
} else if (fc == "LQHPT") {
  args_list[5:9] <- c('linear=TRUE', 'quadratic=TRUE', 'product=TRUE', 'threshold=TRUE', 'hinge=TRUE')
}

# check if it is correct
args_list

# maxent
xm <- maxent(x=predictors, p=pres_train, a=NULL, removeDuplicates=TRUE, nbg=10000,
             factors='biome', path=WD, args=args_list)

# check the results
plot(xm)
response(xm)
e <- evaluate(pres_test, backg, xm, predictors)
e

threshold(e)
plot(e, 'ROC')

# predict the modern case
px <- predict(xm, predictors, progress='text', overwrite=T,
              path=WD, filename=paste0(taxaname,'.grd'))

par(mfrow=c(1,2))
plot(px, main='Maxent, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
plot(px > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')



