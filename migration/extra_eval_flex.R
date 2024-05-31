#extra_eval
training_data = occs
pr_ab = "p"
projection_data = coords_xy[1:100,]
metric = "euclidean" 
univar_comb = FALSE
n_cores = 1
aggreg_factor = 1 

function (training_data, pr_ab, projection_data, metric = "mahalanobis", 
          univar_comb = FALSE, n_cores = 1, aggreg_factor = 1) 
{
  Value <- val <- . <- x <- extrapolation <- NULL
  if (!metric %in% c("euclidean", "mahalanobis")) {
    stop("metric argument must be used with 'euclidean' or 'mahalanobis'")
  }
  if (length(metric) > 1) {
    stop("metric argument must be used with 'euclidean or 'mahalanobis'")
  }
  if (any("data.frame" == class(training_data))) {
    training_data_pr_ab <- training_data[c(names(projection_data), 
                                           pr_ab)] %>% na.omit()
    training_data <- training_data_pr_ab[names(projection_data)] %>% 
      na.omit()
  }
  v0 <- unique(c(names(training_data), names(projection_data)))
  v0 <- sort(v0)
  if (!all(c(all(names(training_data) %in% names(projection_data)), 
             all(names(projection_data) %in% names(training_data))))) {
    stop("colnames of dataframes of env_records, training_data, and projection_data\n        do not match each other", 
         "\nraster layers names:", "\n", paste(sort(unique(unlist(v0))), 
                                               collapse = "\n"))
  }
  if (any("data.frame" %in% class(training_data))) {
    training_data <- training_data[v0]
  }
  else {
    training_data <- training_data[[v0]]
  }
  if (any("SpatRaster" == class(projection_data))) {
    projection_data <- projection_data[[v0]]
    extraraster <- projection_data[[1]]
    extraraster[!is.na(extraraster)] <- 0
    if (aggreg_factor == 1) {
      aggreg_factor <- NULL
    }
    if (!is.null(aggreg_factor)) {
      disag <- extraraster
      if (any("SpatRaster" == class(training_data))) {
        training_data <- terra::aggregate(training_data, 
                                          fact = aggreg_factor, na.rm = TRUE)
      }
      projection_data <- terra::aggregate(projection_data, 
                                          fact = aggreg_factor, na.rm = TRUE)
      extraraster <- terra::aggregate(extraraster, fact = aggreg_factor, 
                                      na.rm = TRUE)
    }
  }
  else {
    projection_data <- projection_data[v0]
  }
  env_calib2 <- training_data <- terra::as.data.frame(training_data, xy = FALSE, 
                                     na.rm = TRUE)
  env_proj2 <- terra::as.data.frame(projection_data, xy = FALSE, 
                                    na.rm = TRUE)
  if (any("SpatRaster" == class(projection_data))) {
    ncell <- rownames(env_proj2) %>% as.numeric()
  }
  
  s_center <- colMeans(env_calib2)
  s_scale <- apply(env_calib2, 2, stats::sd)
  for (i in 1:ncol(env_calib2)) {
    env_calib2[i] <- (env_calib2[i] - s_center[i])/s_scale[i]
  }
  for (i in 1:ncol(env_proj2)) {
    env_proj2[i] <- (env_proj2[i] - s_center[i])/s_scale[i]
  }
  set <- c(seq(1, nrow(env_proj2), 200), nrow(env_proj2) + 
             1)
  if (n_cores == 1) {
    extra <- lapply(seq_len((length(set) - 1)), function(x) {
      rowset <- set[x]:(set[x + 1] - 1)
      if (metric == "euclidean") {
        envdist <- euc_dist(env_proj2[rowset, v0], env_calib2[v0])
        envdist <- sapply(data.frame(t(envdist)), min)
      }
      if (metric == "mahalanobis") {
        envdist <- mah_dist(x = env_proj2[rowset, v0], 
                            y = env_calib2[v0], cov = stats::cov(env_calib2))
        envdist <- sapply(data.frame(t(envdist)), min)
      }
      return(envdist)
    })
  }
  else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    extra <- foreach::foreach(x = seq_len((length(set) - 
                                             1)), .export = c("euc_dist"), .combine = "c") %dopar% 
      {
        rowset <- set[x]:(set[x + 1] - 1)
        if (metric == "euclidean") {
          envdist <- euc_dist(env_proj2[rowset, v0], 
                              env_calib2[v0])
          envdist <- sapply(data.frame(t(envdist)), 
                            min)
        }
        if (metric == "mahalanobis") {
          envdist <- mah_dist(x = env_proj2[rowset, 
                                            v0], y = env_calib2[v0], cov = stats::cov(env_calib2))
          envdist <- sapply(data.frame(t(envdist)), 
                            min)
        }
        envdist
      }
    parallel::stopCluster(cl)
  }
  extra <- unlist(extra)
  if (any("SpatRaster" == class(projection_data))) {
    env_proj2 <- data.frame(distance = extra)
  }
  else {
    env_proj2 <- data.frame(distance = extra, env_proj2)
  }
  rm(extra)
  if (metric == "euclidean") {
    base_stand_distance <- env_calib2 %>% dplyr::summarise_all(., 
                                                               mean) %>% euc_dist(env_calib2, .) %>% mean()
  }
  if (metric == "mahalanobis") {
    base_stand_distance <- env_calib2 %>% dplyr::summarise_all(., 
                                                               mean) %>% mah_dist(x = env_calib2, y = ., cov = stats::cov(env_calib2)) %>% 
      mean()
  }
  env_proj2 <- data.frame(extrapolation = env_proj2$distance/base_stand_distance * 
                            100, env_proj2) %>% dplyr::select(-distance)
  if (any("SpatRaster" == class(projection_data))) {
    extraraster[ncell] <- env_proj2[, "extrapolation"]
    if (!is.null(aggreg_factor)) {
      extraraster <- terra::resample(x = extraraster, 
                                     y = disag)
      extraraster <- terra::mask(extraraster, disag)
    }
    names(extraraster) <- "extrapolation"
    if (univar_comb) {
      rng <- apply(training_data, 2, range, na.rm = TRUE)
      univar_ext <- projection_data
      for (i in 1:terra::nlyr(projection_data)) {
        univar_ext[[i]] <- (projection_data[v0[i]] <= 
                              rng[, v0[i]][1] | projection_data[v0[i]] >= 
                              rng[, v0[i]][2])
      }
      univar_comb_r <- sum(univar_ext)
      univar_comb_r[univar_comb_r > 0] <- 2
      univar_comb_r[univar_comb_r == 0] <- 1
      names(univar_comb_r) <- "uni_comb"
      extraraster <- c(extraraster, univar_comb_r)
    }
  }
  else {
    for (i in names(s_center)) {
      env_proj2[i] <- env_proj2[i] * s_scale[i] + s_center[i]
    }
    if (univar_comb) {
      rng <- apply(training_data, 2, range, na.rm = TRUE)
      univar_ext <- projection_data
      for (i in 1:ncol(projection_data)) {
        univar_ext[, v0[i]] <- (projection_data[, v0[i]] <= 
                                  rng[, v0[i]][1] | projection_data[, v0[i]] >= 
                                  rng[, v0[i]][2])
      }
      univar_comb_r <- apply(univar_ext, 1, sum)
      univar_comb_r[univar_comb_r > 0] <- 2
      univar_comb_r[univar_comb_r == 0] <- 1
      env_proj2 <- env_proj2 %>% dplyr::mutate(univar_comb_r) %>% 
        dplyr::relocate(extrapolation, univar_comb = univar_comb_r)
    }
    return(dplyr::as_tibble(env_proj2))
  }
  return(extraraster)
}

