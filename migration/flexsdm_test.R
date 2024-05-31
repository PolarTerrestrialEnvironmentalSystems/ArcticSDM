#flexsdm

#https://rdrr.io/github/sjevelazco/flexsdm/f/vignettes/v06_Extrapolation_example.Rmd
#Velazco, S. J. E., Brooke, M. R., De Marco Jr., P., Regan, H. M., & Franklin, J. (2023). How far can I extrapolate my species distribution model? Exploring Shape, a novel method. Ecography, 11, e06992. https://doi.org/10.1111/ecog.06992

rm(list=ls())
#### data ####
library(flexsdm)
library(terra)
library(dplyr)
library(patchwork)

# environmental data
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)
names(somevar) <- c("cwd", "tmn", "aet", "ppt_jja")

# species occurence data (presence-only)
data(hespero)
hespero <- hespero %>% dplyr::select(-id)

# California ecoregions
regions <- system.file("external/regions.tif", package = "flexsdm")
regions <- terra::rast(regions)
regions <- terra::as.polygons(regions)
sp_region <- terra::subset(regions, regions$category == "SCR") # ecoregion where *Hesperocyparis stephensonii* is found

# visualize the species occurrences
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  main = "Hesperocyparis stephensonii occurrences"
)
points(hespero[, c("x", "y")], col = "black", pch = 16)
cols <- rep("gray80", 8)
cols[regions$category == "SCR"] <- "yellow"
terra::inset(
  regions,
  loc = "bottomleft",
  scale = .3,
  col = cols
)

#### Delimit calibration area ####
ca <- calib_area(
  data = hespero,
  x = "x",
  y = "y",
  method = c("buffer", width = 25000),
  crs = crs(somevar)
)

# visualize the species occurrences & calibration area
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  main = "Calibration area and occurrences"
)
plot(ca, add = TRUE)
points(hespero[, c("x", "y")], col = "black", pch = 16)

#### Create pseudo-absence data ####
# Sample the same number of species presences
set.seed(10)
psa <- sample_pseudoabs(
  data = hespero,
  x = "x",
  y = "y",
  n = sum(hespero$pr_ab), # number of pseudo-absence points equal to number of presences
  method = "random",
  rlayer = somevar,
  calibarea = ca
)

# Visualize species presences and pseudo-absences
plot(
  sp_region,
  col = "gray80",
  legend = FALSE,
  axes = FALSE,
  xlim = c(289347, 353284),
  ylim = c(-598052, -520709),
  main = "Presence = yellow, Pseudo-absence = black"
)
plot(ca, add = TRUE)
points(psa[, c("x", "y")], cex = 0.8, pch = 16, col = "black") # Pseudo-absences
points(hespero[, c("x", "y")], col = "yellow", pch = 16, cex = 1.5) # Presences


# Bind a presences and pseudo-absences
hespero_pa <- bind_rows(hespero, psa)
hespero_pa # Presence-Pseudo-absence database

#### Partition data for evaluating models ####
set.seed(10)

# Repeated K-fold method
hespero_pa2 <- part_random(
  data = hespero_pa,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 5, replicates = 10)
)

#### Extracting environmental values ####
hespero_pa3 <-
  sdm_extract(
    data = hespero_pa2,
    x = "x",
    y = "y",
    env_layer = somevar,
    variables = c("cwd", "tmn", "aet", "ppt_jja")
  )

#### Modeling ####
mglm <-
  fit_glm(
    data = hespero_pa3,
    response = "pr_ab",
    predictors = c("cwd", "tmn", "aet", "ppt_jja"),
    partition = ".part",
    thr = "max_sens_spec"
  )

mgbm <- fit_gbm(
  data = hespero_pa3,
  response = "pr_ab",
  predictors = c("cwd", "tmn", "aet", "ppt_jja"),
  partition = ".part",
  thr = "max_sens_spec"
)

msvm <- fit_svm(
  data = hespero_pa3,
  response = "pr_ab",
  predictors = c("cwd", "tmn", "aet", "ppt_jja"),
  partition = ".part",
  thr = "max_sens_spec"
)


mpred <- sdm_predict(
  models = list(mglm, mgbm, msvm),
  pred = somevar,
  con_thr = TRUE,
  predict_area = NULL
)

#### Comparing our models ####
par(mfrow = c(1, 3))
plot(mpred$glm, main = "GLM")
# points(hespero$x, hespero$y, pch = 19)
plot(mpred$gbm, main = "GBM")
# points(hespero$x, hespero$y, pch = 19)
plot(mpred$svm, main = "SVM")
# points(hespero$x, hespero$y, pch = 19)

#### Partial dependence plots to explore the impact of predictor conditions on suitability ####

#Uni and bivariate partial dependence plots for the GLM:
p_pdp(model = mglm$model, training_data = hespero_pa3, projection_data = somevar)
p_bpdp(model = mglm$model, training_data = hespero_pa3, training_boundaries = "convexh")

#Uni and bivariate partial dependence plots for the GBM:
p_pdp(model = mgbm$model, training_data = hespero_pa3, projection_data = somevar)
p_bpdp(model = mgbm$model, training_data = hespero_pa3, training_boundaries = "convexh", resolution = 100)

#Uni and bivariate partial dependence plots for the SVM:
p_pdp(model = msvm$model, training_data = hespero_pa3, projection_data = somevar)
p_bpdp(model = msvm$model, training_data = hespero_pa3, training_boundaries = "convexh")

#### Extrapolation evaluation ####
#Using Mahalanobis distance:
  
  xp_m <-
  extra_eval(
    training_data = hespero_pa3,
    pr_ab = "pr_ab",
    projection_data = somevar,
    metric = "mahalanobis",
    univar_comb = TRUE,
    n_cores = 1,
    aggreg_factor = 1
  )
xp_m

#The output of the extra_eval function is a SpatRaster, 
#showing the degree of extrapolation across the projection area, .
#as estimated by the Shape method.

cl <- c("#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")

par(mfrow = c(1, 2))
plot(xp_m$extrapolation, main = "Shape metric", col = cl)
plot(xp_m$uni_comb, main = "Univariate (1) and \n combinatorial (2) extrapolation", col = cl)

#We can also explore extrapolation or suitability patterns in environmental and geographic space, 
#using just one function. To do that, we will use the p_extra function. This function plots a ggplot object.

#Let's start with our extrapolation evaluation. 
#These plots show that areas with high extrapolation (dark blue) are far from the training data (shown in black) 
#in both environmental and geographic space.

#The higher extrapolation values extrapolation area in the northwestern portion of the CFP.

p_extra(
  training_data = hespero_pa3,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  color_p = "black",
  extra_suit_data = xp_m,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05
)

#Let's explore univariate and combinatorial extrapolation. 
#The former is defined as the projecting data outside range of training conditions, 
#while the combinatorial extrapolation area those projecting data within the range of training conditions.

p_extra(
  training_data = hespero_pa3,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  color_p = "black",
  extra_suit_data = xp_m$uni_comb,
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
  color_gradient = c("#B3DC2B", "#30678D"),
  alpha_p = 0.2
)

#### Truncating SDMs predictions based on extrapolation thresholds ####
p_extra(
  training_data = hespero_pa3,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  color_p = "black",
  extra_suit_data = as.numeric(xp_m$extrapolation < 50),
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
  color_gradient = c("gray", "#FDE725"),
  alpha_p = 0.5
) + plot_annotation(subtitle = "Binary extrapolation pattern with using a threshold of 50")

p_extra(
  training_data = hespero_pa3,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  color_p = "black",
  extra_suit_data = as.numeric(xp_m$extrapolation < 100),
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
  color_gradient = c("gray", "#FDE725"),
  alpha_p = 0.5
) + plot_annotation(subtitle = "Binary extrapolation pattern with using a threshold of 100")

p_extra(
  training_data = hespero_pa3,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  color_p = "black",
  extra_suit_data = as.numeric(xp_m$extrapolation < 500),
  projection_data = somevar,
  geo_space = TRUE,
  prop_points = 0.05,
  color_gradient = c("gray", "#FDE725"),
  alpha_p = 0.5
) + plot_annotation(subtitle = "Binary extrapolation pattern with using a threshold of 500")

#### truncating ####
glm_trunc <- extra_truncate(
  suit = mpred$glm,
  extra = xp_m,
  threshold = c(50, 100, 500),
  trunc_value = 0
)

gbm_trunc <- extra_truncate(
  suit = mpred$gbm,
  extra = xp_m,
  threshold = c(50, 100, 500),
  trunc_value = 0
)

svm_trunc <- extra_truncate(
  suit = mpred$svm,
  extra = xp_m,
  threshold = c(50, 100, 500),
  trunc_value = 0
)

par(mfrow = c(3, 3))
plot(glm_trunc$`50`, main = "GLM; extra threshold = 50", col = cl)
plot(glm_trunc$`100`, main = "GLM; extra threshold = 100", col = cl)
plot(glm_trunc$`500`, main = "GLM; extra threshold = 500", col = cl)
plot(gbm_trunc$`50`, main = "GBM; extra threshold = 50", col = cl)
plot(gbm_trunc$`100`, main = "GBM; extra threshold = 100", col = cl)
plot(gbm_trunc$`500`, main = "GBM; extra threshold = 500", col = cl)
plot(svm_trunc$`50`, main = "SVM; extra threshold = 50", col = cl)
plot(svm_trunc$`100`, main = "SVM; extra threshold = 100", col = cl)
plot(svm_trunc$`500`, main = "SVM; extra threshold = 500", col = cl)

