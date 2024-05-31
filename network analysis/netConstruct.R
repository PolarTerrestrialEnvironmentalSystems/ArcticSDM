#netConstruct

data = netwide0
verbose=2

library(NetCoMi)

netConstruct <- function (data, data2 = NULL, dataType = "counts", group = NULL, 
          matchDesign = NULL, taxRank = NULL, measure = "spieceasi", 
          measurePar = NULL, jointPrepro = NULL, filtTax = "none", 
          filtTaxPar = NULL, filtSamp = "none", filtSampPar = NULL, 
          zeroMethod = "none", zeroPar = NULL, normMethod = "none", 
          normPar = NULL, sparsMethod = "t-test", thresh = 0.3, alpha = 0.05, 
          adjust = "adaptBH", trueNullMethod = "convest", lfdrThresh = 0.2, 
          nboot = 1000L, assoBoot = NULL, cores = 1L, logFile = "log.txt", 
          softThreshType = "signed", softThreshPower = NULL, softThreshCut = 0.8, 
          kNeighbor = 3L, knnMutual = FALSE, dissFunc = "signed", dissFuncPar = NULL, 
          simFunc = NULL, simFuncPar = NULL, scaleDiss = TRUE, weighted = TRUE, 
          sampleSize = NULL, verbose = 2, seed = NULL) {
  argsIn <- as.list(environment())
  assoType <- distNet <- needfrac <- needint <- NULL
  if (verbose %in% 2:3) {
    message("Checking input arguments ... ", appendLF = FALSE)
  }
  argsOut <- .checkArgsNetConst(argsIn)
  for (i in 1:length(argsOut)) {
    assign(names(argsOut)[i], argsOut[[i]])
  }
  if (verbose %in% 2:3) 
    message("Done.")
  if (dataType == "phyloseq") {
    dataType <- "counts"
  }
  if (dataType == "counts") {
    if (inherits(data, "phyloseq")) {
      otutab <- phyloseq::otu_table(data)
      if (!is.null(taxRank)) {
        taxtab <- as(tax_table(data), "matrix")
      }
      if (attributes(otutab)$taxa_are_rows) {
        data <- t(as(otutab, "matrix"))
      }
      else {
        data <- as(otutab, "matrix")
      }
      if (!is.null(taxRank)) {
        if (!taxRank %in% colnames(taxtab)) {
          stop("Argument \"taxRank\" must match column names of the taxonomic ", 
               "table:\n", paste(colnames(taxtab), collapse = ", "))
        }
        if (any(duplicated(taxtab[, taxRank]))) {
          stop(paste0("Taxa names of chosen taxonomic rank must be unique. ", 
                      "Consider using NetCoMi's function renameTaxa()."))
        }
        colnames(data) <- taxtab[, taxRank]
      }
    }
    else {
      if (is.null(rownames(data))) {
        message("Row names are numbered because sample names were missing.\n")
        rownames(data) <- 1:nrow(data)
      }
      if (is.null(colnames(data))) {
        message("Column names are numbered because taxa names were missing.\n")
        colnames(data) <- 1:ncol(data)
      }
      if (identical(colnames(data), rownames(data))) {
        warning(paste0("Row names and column names of 'data' are equal. ", 
                       "Ensure 'data' is a count matrix."))
      }
    }
    if (!is.null(data2)) {
      if (inherits(data2, "phyloseq")) {
        otutab <- phyloseq::otu_table(data2)
        if (!is.null(taxRank)) {
          taxtab <- as(tax_table(data2), "matrix")
        }
        if (attributes(otutab)$taxa_are_rows) {
          data2 <- t(as(otutab, "matrix"))
        }
        else {
          data2 <- as(otutab, "matrix")
        }
        if (!is.null(taxRank)) {
          if (!taxRank %in% colnames(taxtab)) {
            stop("Argument \"taxRank\" must match column names of the taxonomic ", 
                 "table:\n", paste(colnames(taxtab), collapse = ", "))
          }
          if (any(duplicated(taxtab[, taxRank]))) {
            stop(paste0("Taxa names of chosen taxonomic rank must be unique. ", 
                        "Consider using NetCoMi's function renameTaxa()."))
          }
          colnames(data2) <- taxtab[, taxRank]
        }
      }
      if (is.null(rownames(data2))) {
        message("Row names of 'data2' are numbered because sample names ", 
                "were missing.\n")
        rownames(data2) <- 1:nrow(data2)
      }
      if (is.null(colnames(data2))) {
        message("Column names of 'data2' are numbered because taxa names ", 
                "were missing.\n")
        colnames(data2) <- 1:ncol(data2)
      }
      if (identical(colnames(data2), rownames(data2))) {
        warning(paste0("Row names and column names of 'data2' are equal. ", 
                       "Ensure 'data2' is a count matrix."))
      }
    }
  }
  else {
    if (xor(is.null(rownames(data)), is.null(colnames(data))) || 
        !identical(rownames(data), colnames(data))) {
      stop("Row and column names must match.")
    }
    if (is.null(rownames(data))) {
      message("Row and columns names of 'data' missing. ", 
              "Numbers are used instead.")
      rownames(data) <- colnames(data) <- 1:nrow(data)
    }
    if (!is.null(data2)) {
      if (xor(is.null(rownames(data2)), is.null(colnames(data2))) || 
          !identical(rownames(data2), colnames(data2))) {
        stop("Row and column names must match.")
      }
      if (is.null(rownames(data2))) {
        message("Row and columns names of 'data2' missing. ", 
                "Numbers are used instead.")
        rownames(data2) <- colnames(data2) <- 1:nrow(data2)
      }
    }
  }
  plausCheck <- .checkPlausNetConst(dataType = dataType, assoType = assoType, 
                                    data2 = data2, measure = measure, normMethod = normMethod, 
                                    zeroMethod = zeroMethod, sparsMethod = sparsMethod, dissFunc = dissFunc, 
                                    sampleSize = sampleSize, verbose = verbose)
  for (i in 1:length(plausCheck)) {
    assign(names(plausCheck)[i], plausCheck[[i]])
  }
  .checkPackNetConst(measure = measure, zeroMethod = zeroMethod, 
                     normMethod = normMethod, sparsMethod = sparsMethod, adjust = adjust)
  if (!is.null(seed)) 
    set.seed(seed)
  twoNets <- ifelse(is.null(data2) & is.null(group), FALSE, 
                    TRUE)
  if (twoNets) {
    if (!is.null(group) && !is.null(data2)) {
      stop("Only one of the arguments 'group' and 'data2' may be defined.")
    }
    if (!is.null(group)) {
      if (is.null(jointPrepro)) {
        if (distNet) {
          jointPrepro <- FALSE
        }
        else {
          jointPrepro <- TRUE
        }
      }
    }
    else if (!is.null(data2) && is.null(jointPrepro)) {
      jointPrepro <- FALSE
    }
    if (jointPrepro && distNet) {
      stop("'jointPrepro' is TRUE but data must be normalized separately ", 
           "for dissimilarity measures.")
    }
    if (length(thresh) == 1) 
      thresh <- c(thresh, thresh)
    if (length(alpha) == 1) 
      alpha <- c(alpha, alpha)
    if (length(lfdrThresh) == 1) 
      lfdrThresh <- c(lfdrThresh, lfdrThresh)
    if (length(softThreshPower) == 1) 
      softThreshPower <- c(softThreshPower, softThreshPower)
    if (length(softThreshCut) == 1) 
      softThreshCut <- c(softThreshCut, softThreshCut)
    if (length(sampleSize) == 1) 
      sampleSize <- c(sampleSize, sampleSize)
  }
  if (dataType == "counts") {
    if (!(filtTax[1] == "none" & filtSamp[1] == "none") & 
        verbose %in% 1:3) {
      message("Data filtering ...")
    }
    countMat1 <- t(apply(data, 1, function(x) as.numeric(x)))
    colnames(countMat1) <- colnames(data)
    if (!is.null(data2)) {
      countMat2 <- t(apply(data2, 1, function(x) as.numeric(x)))
      colnames(countMat2) <- colnames(data2)
    }
    countMatJoint <- countsJointOrig <- NULL
    if (twoNets) {
      if (!is.null(group)) {
        if (distNet) {
          if (!(is.vector(group) || is.factor(group))) {
            stop("'group' must be of type vector or factor.")
          }
          if (length(group) != ncol(countMat1)) {
            stop("Length of 'group' must match the number of columns of 'data'.")
          }
          group <- as.numeric(group)
          if (is.null(names(group))) {
            names(group) <- colnames(countMat1)
          }
          else {
            if (!all(colnames(countMat1) %in% names(group))) {
              stop("Names of 'group' must match column names of 'data'.")
            }
          }
          if (any(is.na(countMat1))) {
            countMat1 <- countMat1[complete.cases(countMat1), 
                                   , drop = FALSE]
            if (verbose %in% 1:3) 
              message("Samples with NAs removed.")
          }
        }
        else {
          if (!(is.vector(group) || is.factor(group))) {
            stop("'group' must be of type vector or factor.")
          }
          if (length(group) != nrow(countMat1)) {
            stop("Length of 'group' must match the number of rows of 'data'.")
          }
          if (is.character(group)) {
            group <- as.factor(group)
          }
          group <- as.numeric(group)
          if (is.null(names(group))) {
            names(group) <- rownames(countMat1)
          }
          else {
            if (!all(rownames(countMat1) %in% names(group))) {
              stop("Names of 'group' must match column names of 'data'.")
            }
          }
          if (any(is.na(countMat1))) {
            if (!is.null(matchDesign)) {
              stop("Data set contains NAs. ", "Cannot be removed if a matched-group design is used.")
            }
            data_tmp <- cbind(countMat1, group)
            data_tmp <- data_tmp[complete.cases(data_tmp), 
                                 , drop = FALSE]
            countMat1 <- data_tmp[, 1:(ncol(countMat1))]
            group <- data_tmp[, ncol(data_tmp)]
            if (verbose %in% 1:3) 
              message("Samples with NAs removed.")
          }
          if (!is.null(matchDesign)) {
            if (!(matchDesign[2]/matchDesign[1]) == (table(group)[2]/table(group)[1])) {
              stop("Group vector not consistent with matched-group design.")
            }
          }
        }
        if (jointPrepro) {
          countMatJoint <- countMat1
          countMat2 <- NULL
        }
        else {
          if (distNet) {
            splitcount <- split(as.data.frame(t(countMat1)), 
                                as.factor(group))
            if (length(splitcount) != 2) 
              stop("Argument 'group' has to be binary.")
            groups <- names(splitcount)
            countMat1 <- t(as.matrix(splitcount[[1]]))
            countMat2 <- t(as.matrix(splitcount[[2]]))
          }
          else {
            splitcount <- split(as.data.frame(countMat1), 
                                as.factor(group))
            if (length(splitcount) != 2) 
              stop("Argument 'group' has to be binary.")
            groups <- names(splitcount)
            countMat1 <- as.matrix(splitcount[[1]])
            countMat2 <- as.matrix(splitcount[[2]])
          }
        }
      }
      else {
        groups <- c("1", "2")
        if (distNet) {
          if (!identical(colnames(countMat1), colnames(countMat2))) {
            if (!all(rownames(countMat1) %in% rownames(countMat2))) {
              if (verbose > 0) 
                message("Intersection of samples selected.")
            }
            sel <- intersect(rownames(countMat1), rownames(countMat2))
            if (length(sel) == 0) 
              stop("Data sets contain different samples")
            countMat1 <- countMat1[sel, , drop = FALSE]
            countMat2 <- countMat2[sel, , drop = FALSE]
          }
          if (any(is.na(countMat1)) || any(is.na(countMat2))) {
            if (verbose %in% 1:3) 
              message("Samples with NAs removed.")
            keep <- intersect(which(complete.cases(countMat1)), 
                              which(complete.cases(countMat2)))
            countMat1 <- countMat1[keep, , drop = FALSE]
            countMat2 <- countMat2[keep, , drop = FALSE]
          }
        }
        else {
          if (!identical(colnames(countMat1), colnames(countMat2))) {
            if (!all(colnames(countMat1) %in% colnames(countMat2))) {
              if (verbose > 0) 
                message("Intersection of taxa selected.")
            }
            sel <- intersect(colnames(countMat1), colnames(countMat2))
            if (length(sel) == 0) 
              stop("Data sets contain different taxa.")
            countMat1 <- countMat1[, sel, drop = FALSE]
            countMat2 <- countMat2[, sel, drop = FALSE]
          }
          if (any(is.na(countMat1)) || any(is.na(countMat2))) {
            if (!is.null(matchDesign)) {
              stop("Data contain NAs. ", "Cannot be removed if a matched-group design is used.")
            }
            if (verbose %in% 1:3) 
              message("Samples with NAs removed.")
            countMat1 <- countMat1[complete.cases(countMat1), 
                                   , drop = FALSE]
            countMat2 <- countMat2[complete.cases(countMat2), 
                                   , drop = FALSE]
          }
          if (!is.null(matchDesign)) {
            if (!(matchDesign[2]/matchDesign[1]) == (nrow(countMat2)/nrow(countMat1))) {
              stop("Sample sizes not consistent with matched-group design.")
            }
          }
        }
        if (jointPrepro) {
          if (any(rownames(countMat2) %in% rownames(countMat1))) {
            if (verbose %in% 1:3) {
              message("\"*\" Added to duplicated sample names in group 2.")
            }
            dupidx <- which(rownames(countMat2) %in% 
                              rownames(countMat1))
            rownames(countMat2)[dupidx] <- paste0(rownames(countMat2)[dupidx], 
                                                  "*")
          }
          countMatJoint <- rbind(countMat1, countMat2)
          n1 <- nrow(countMat1)
          n2 <- nrow(countMat2)
        }
      }
    }
    else {
      if (any(is.na(data))) {
        if (verbose %in% 1:3) 
          message("Samples with NAs removed.")
        data <- data[complete.cases(data), , drop = FALSE]
      }
      countMatJoint <- countMat1
      countMat2 <- NULL
    }
    if (!twoNets || jointPrepro) {
      keepRows <- .filterSamples(countMat = countMatJoint, 
                                 filter = filtSamp, filterParam = filtSampPar)
      if (length(keepRows) == 0) {
        stop("No samples remaining after filtering.")
      }
      if (length(keepRows) != nrow(countMatJoint)) {
        n_old <- nrow(countMatJoint)
        countMatJoint <- countMatJoint[keepRows, , drop = FALSE]
        if (!distNet) 
          group <- group[keepRows]
        if (verbose %in% 2:3) {
          message(n_old - nrow(countMatJoint), " samples removed.")
        }
      }
    }
    else {
      n_old1 <- nrow(countMat1)
      n_old2 <- nrow(countMat2)
      keepRows1 <- .filterSamples(countMat = countMat1, 
                                  filter = filtSamp, filterParam = filtSampPar)
      keepRows2 <- .filterSamples(countMat = countMat2, 
                                  filter = filtSamp, filterParam = filtSampPar)
      if (distNet) {
        keepRows <- intersect(keepRows1, keepRows2)
        countMat1 <- countMat1[keepRows, , drop = FALSE]
        countMat2 <- countMat2[keepRows, , drop = FALSE]
        if (n_old1 - nrow(countMat1) != 0 && verbose %in% 
            2:3) {
          message(n_old1 - nrow(countMat1), " samples removed in each data set.")
        }
        if (length(keepRows) == 0) {
          stop("No samples remaining after filtering and building the ", 
               "intercept of the remaining samples.")
        }
      }
      else {
        countMat1 <- countMat1[keepRows1, , drop = FALSE]
        countMat2 <- countMat2[keepRows2, , drop = FALSE]
        if (n_old1 - nrow(countMat1) != 0 || n_old2 - 
            nrow(countMat2)) {
          if (verbose %in% 2:3) {
            message(n_old1 - nrow(countMat1), " samples removed in data set 1.")
            message(n_old2 - nrow(countMat2), " samples removed in data set 2.")
          }
        }
      }
      if (length(keepRows1) == 0) {
        stop("No samples remaining in group 1 after filtering.")
      }
      if (length(keepRows2) == 0) {
        stop("No samples remaining in group 2 after filtering.")
      }
    }
    if (!twoNets || jointPrepro) {
      keepCols <- .filterTaxa(countMat = countMatJoint, 
                              filter = filtTax, filterParam = filtTaxPar)
      if (length(keepCols) != ncol(countMatJoint)) {
        p_old <- ncol(countMatJoint)
        countMatJoint <- countMatJoint[, keepCols, drop = FALSE]
        if (verbose %in% 2:3) 
          message(p_old - ncol(countMatJoint), " taxa removed.")
        if (distNet) 
          group <- group[keepCols]
      }
      if (length(keepCols) == 0) {
        stop("No taxa remaining after filtering.")
      }
    }
    else {
      p_old1 <- ncol(countMat1)
      p_old2 <- ncol(countMat2)
      keepCols1 <- .filterTaxa(countMat = countMat1, filter = filtTax, 
                               filterParam = filtTaxPar)
      keepCols2 <- .filterTaxa(countMat = countMat2, filter = filtTax, 
                               filterParam = filtTaxPar)
      if (!distNet) {
        keepCols <- intersect(keepCols1, keepCols2)
        countMat1 <- countMat1[, keepCols, drop = FALSE]
        countMat2 <- countMat2[, keepCols, drop = FALSE]
        if (length(keepCols) == 0) {
          stop("No taxa remaining after filtering and building the intercept ", 
               "of the remaining taxa.")
        }
        if (p_old1 - dim(countMat1)[2] != 0 && verbose %in% 
            2:3) {
          message(p_old1 - dim(countMat1)[2], " taxa removed in each data set.")
        }
      }
      else {
        countMat1 <- countMat1[, keepCols1, drop = FALSE]
        countMat2 <- countMat2[, keepCols2, drop = FALSE]
        if (p_old1 - dim(countMat1)[2] != 0 || p_old2 - 
            dim(countMat2)[2] != 0) {
          if (verbose %in% 2:3) {
            message(p_old1 - dim(countMat1)[2], " taxa removed in data set 1.")
            message(p_old2 - dim(countMat2)[2], " taxa removed in data set 2.")
          }
        }
      }
      if (length(keepCols1) == 0) {
        stop("No samples remaining in group 1 after filtering.")
      }
      if (length(keepCols2) == 0) {
        stop("No samples remaining in group 2 after filtering.")
      }
    }
    rmZeroSum <- function(countMat, group, matchDesign) {
      rs <- Matrix::rowSums(countMat)
      if (any(rs == 0)) {
        rmRows <- which(rs == 0)
        if (!is.null(matchDesign)) {
          stop(paste0("The following samples have an overall sum of zero ", 
                      "but cannot be removed if a matched-group design is used: "), 
               paste0(rmRows, sep = ","))
        }
        countMat <- countMat[-rmRows, , drop = FALSE]
        if (!is.null(group) & !distNet & length(rmRows) != 
            0) {
          group <- group[-rmRows]
        }
      }
      return(list(countMat = countMat, group = group))
    }
    if (!twoNets || jointPrepro) {
      n_old <- nrow(countMatJoint)
      rmZeroSum_res <- rmZeroSum(countMatJoint, group = group, 
                                 matchDesign = matchDesign)
      countMatJoint <- rmZeroSum_res$countMat
      group <- rmZeroSum_res$group
      if (verbose %in% 2:3 && n_old != nrow(countMatJoint)) {
        message(paste(n_old - nrow(countMatJoint), "rows with zero sum removed."))
      }
    }
    else {
      n_old1 <- nrow(countMat1)
      n_old2 <- nrow(countMat2)
      countMat1 <- rmZeroSum(countMat1, group = group, 
                             matchDesign = matchDesign)$countMat
      countMat2 <- rmZeroSum(countMat2, group = group, 
                             matchDesign = matchDesign)$countMat
      if (distNet) {
        keep <- intersect(rownames(countMat1), rownames(countMat2))
        countMat1 <- countMat1[keep, , drop = FALSE]
        countMat2 <- countMat2[keep, , drop = FALSE]
        if (verbose %in% 2:3 && n_old1 != nrow(countMat1)) {
          message(paste(n_old1 - nrow(countMat1), "rows with zero sum removed in both groups."))
        }
      }
      else {
        if (verbose %in% 2:3) {
          if (n_old1 != nrow(countMat1)) {
            message(paste(n_old1 - nrow(countMat1), "rows with zero sum removed in group 1."))
          }
          if (n_old2 != nrow(countMat2)) {
            message(paste(n_old2 - nrow(countMat2), "rows with zero sum removed in group 2."))
          }
        }
      }
    }
    if (!twoNets || jointPrepro) {
      if (verbose %in% 1:3) 
        message(ncol(countMatJoint), " taxa and ", nrow(countMatJoint), 
                " samples remaining.")
    }
    else {
      if (verbose %in% 1:3) {
        message(ncol(countMat1), " taxa and ", nrow(countMat1), 
                " samples remaining in group 1.")
        message(ncol(countMat2), " taxa and ", nrow(countMat2), 
                " samples remaining in group 2.")
      }
    }
    if (twoNets) {
      if (jointPrepro) {
        countsJointOrig <- countMatJoint
        attributes(countsJointOrig)$scale <- "counts"
        if (!is.null(data2)) {
          n1 <- sum(rownames(countMatJoint) %in% rownames(countMat1))
          n2 <- sum(rownames(countMatJoint) %in% rownames(countMat2))
        }
        countsOrig1 <- countsOrig2 <- NULL
      }
      else {
        countsOrig1 <- countMat1
        countsOrig2 <- countMat2
        attributes(countsOrig1)$scale <- "counts"
        attributes(countsOrig2)$scale <- "counts"
      }
    }
    else {
      countsOrig1 <- countMatJoint
      attributes(countsOrig1)$scale <- "counts"
      countsOrig2 <- NULL
    }
    if (zeroMethod != "none") {
      if (!twoNets || jointPrepro) {
        if (verbose %in% 2:3) 
          message("\nZero treatment:")
        countMatJoint <- .zeroTreat(countMat = countMatJoint, 
                                    zeroMethod = zeroMethod, zeroParam = zeroPar, 
                                    needfrac = needfrac, needint = needint, verbose = verbose)
      }
      else {
        if (verbose %in% 2:3) 
          message("\nZero treatment in group 1:")
        countMat1 <- .zeroTreat(countMat = countMat1, 
                                zeroMethod = zeroMethod, zeroParam = zeroPar, 
                                needfrac = needfrac, needint = needint, verbose = verbose)
        if (verbose %in% 2:3) 
          message("\nZero treatment in group 2:")
        countMat2 <- .zeroTreat(countMat = countMat2, 
                                zeroMethod = zeroMethod, zeroParam = zeroPar, 
                                needfrac = needfrac, needint = needint, verbose = verbose)
      }
    }
    else {
      attributes(countMat1)$scale <- "counts"
      if (!is.null(countMat2)) {
        attributes(countMat2)$scale <- "counts"
      }
      if (!is.null(countMatJoint)) {
        attributes(countMatJoint)$scale <- "counts"
      }
    }
    if (!twoNets || jointPrepro) {
      if (verbose %in% 2:3 & (normMethod != "none" || needfrac)) {
        message("\nNormalization:")
      }
      countMatJoint <- .normCounts(countMat = countMatJoint, 
                                   normMethod = normMethod, normParam = normPar, 
                                   zeroMethod = zeroMethod, needfrac = needfrac, 
                                   verbose = verbose)
      if (!twoNets) {
        counts1 <- countMatJoint
        counts2 <- NULL
        groups <- NULL
        sampleSize <- nrow(counts1)
      }
    }
    else {
      if (verbose %in% 2:3 & (normMethod != "none" || needfrac)) {
        message("\nNormalization in group 1:")
      }
      countMat1 <- .normCounts(countMat = countMat1, normMethod = normMethod, 
                               normParam = normPar, zeroMethod = zeroMethod, 
                               needfrac = needfrac, verbose = verbose)
      if (verbose %in% 2:3 & (normMethod != "none" || needfrac)) {
        message("\nNormalization in group 2:")
      }
      countMat2 <- .normCounts(countMat = countMat2, normMethod = normMethod, 
                               normParam = normPar, zeroMethod = zeroMethod, 
                               needfrac = needfrac, verbose = verbose)
    }
    if (twoNets) {
      if (jointPrepro) {
        if (!is.null(group)) {
          splitcount <- split(as.data.frame(countMatJoint), 
                              as.factor(group))
          if (length(splitcount) != 2) 
            stop("Argument 'group' has to be binary.")
          groups <- names(splitcount)
          counts1 <- as.matrix(splitcount[[1]])
          counts2 <- as.matrix(splitcount[[2]])
          sampleSize <- c(nrow(counts1), nrow(counts2))
          rm(countMatJoint, countMat1)
        }
        else {
          counts1 <- countMatJoint[1:n1, , drop = FALSE]
          counts2 <- countMatJoint[(n1 + 1):(n1 + n2), 
                                   , drop = FALSE]
          sampleSize <- c(nrow(counts1), nrow(counts2))
        }
      }
      else {
        counts1 <- countMat1
        counts2 <- countMat2
        if (distNet) {
          keep <- intersect(rownames(counts1), rownames(counts2))
          counts1 <- counts1[keep, , drop = FALSE]
          counts2 <- counts2[keep, , drop = FALSE]
        }
        else {
          keep <- intersect(colnames(counts1), colnames(counts2))
          counts1 <- counts1[, keep, drop = FALSE]
          counts2 <- counts2[, keep, drop = FALSE]
        }
      }
    }
    else {
      rm(countMat1, countMatJoint)
    }
    if (verbose %in% 2:3) {
      if (distNet) {
        txt.tmp <- "dissimilarities"
      }
      else {
        txt.tmp <- "associations"
      }
      message("\nCalculate '", measure, "' ", txt.tmp, 
              " ... ", appendLF = FALSE)
    }
    assoMat1 <- .calcAssociation(countMat = counts1, measure = measure, 
                                 measurePar = measurePar, verbose = verbose)
    if (verbose %in% 2:3) 
      message("Done.")
    if (twoNets) {
      if (verbose %in% 2:3) {
        message("\nCalculate ", txt.tmp, " in group 2 ... ", 
                appendLF = FALSE)
      }
      assoMat2 <- .calcAssociation(countMat = counts2, 
                                   measure = measure, measurePar = measurePar, verbose = verbose)
      if (verbose %in% 2:3) 
        message("Done.")
    }
    else {
      assoMat2 <- NULL
    }
  }
  else {
    assoMat1 <- data
    assoMat2 <- data2
    counts1 <- NULL
    counts2 <- NULL
    countsJointOrig <- NULL
    countsOrig1 <- NULL
    countsOrig2 <- NULL
    groups <- NULL
  }
  if (distNet) {
    dissEst1 <- assoMat1
    if (any(is.infinite(dissEst1)) & scaleDiss) {
      scaleDiss <- FALSE
      warning("Dissimilarity matrix contains infinite values and cannot be ", 
              "scaled to [0,1].")
    }
    if (scaleDiss) {
      assoUpper <- assoMat1[upper.tri(assoMat1)]
      if (length(assoUpper) == 1) {
        warning("Network consists of only two nodes.")
        assoMat1[1, 2] <- assoMat1[2, 1] <- 1
      }
      else {
        assoMat1 <- (assoMat1 - min(assoUpper))/(max(assoUpper) - 
                                                   min(assoUpper))
        diag(assoMat1) <- 0
      }
    }
    dissScale1 <- assoMat1
    if (verbose %in% 2:3) {
      if (sparsMethod != "none") {
        message("\nSparsify dissimilarities via '", sparsMethod, 
                "' ... ", appendLF = FALSE)
      }
    }
    sparsReslt <- .sparsify(assoMat = assoMat1, countMat = counts1, 
                            sampleSize = sampleSize[1], measure = measure, measurePar = measurePar, 
                            assoType = assoType, sparsMethod = sparsMethod, thresh = thresh[1], 
                            alpha = alpha[1], adjust = adjust, lfdrThresh = lfdrThresh[1], 
                            trueNullMethod = trueNullMethod, nboot = nboot, assoBoot = assoBoot, 
                            softThreshType = softThreshType, softThreshPower = softThreshPower[1], 
                            softThreshCut = softThreshCut[1], cores = cores, 
                            logFile = logFile, kNeighbor = kNeighbor, knnMutual = knnMutual, 
                            verbose = verbose, seed = seed)
    if (verbose %in% 2:3 & sparsMethod != "none") 
      message("Done.")
    assoMat1 <- NULL
    assoEst1 <- NULL
    dissMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    simMat1 <- .transToSim(x = dissMat1, simFunc = simFunc, 
                           simFuncPar = simFuncPar)
    adjaMat1 <- .transToAdja(x = simMat1, weighted = weighted)
    if (twoNets) {
      dissEst2 <- assoMat2
      if (scaleDiss) {
        assoUpper <- assoMat2[upper.tri(assoMat2)]
        if (length(assoUpper) == 1) {
          warning("Network consists of only two nodes.")
          assoMat2[1, 2] <- assoMat2[2, 1] <- 1
        }
        else {
          assoMat2 <- (assoMat2 - min(assoUpper))/(max(assoUpper) - 
                                                     min(assoUpper))
          diag(assoMat2) <- 0
        }
      }
      dissScale2 <- assoMat2
      if (verbose %in% 2:3) {
        if (sparsMethod != "none") {
          message("\nSparsify dissimilarities in group 2 ... ", 
                  appendLF = FALSE)
        }
      }
      sparsReslt <- .sparsify(assoMat = assoMat2, countMat = counts2, 
                              sampleSize = sampleSize[2], measure = measure, 
                              measurePar = measurePar, assoType = assoType, 
                              sparsMethod = sparsMethod, thresh = thresh[2], 
                              alpha = alpha[2], adjust = adjust, lfdrThresh = lfdrThresh[2], 
                              trueNullMethod = trueNullMethod, nboot = nboot, 
                              assoBoot = assoBoot, softThreshType = softThreshType, 
                              softThreshPower = softThreshPower[2], softThreshCut = softThreshCut[2], 
                              cores = cores, logFile = logFile, kNeighbor = kNeighbor, 
                              knnMutual = knnMutual, verbose = verbose, seed = seed)
      if (verbose %in% 2:3 & sparsMethod != "none") 
        message("Done.")
      assoMat2 <- NULL
      assoEst2 <- NULL
      dissMat2 <- sparsReslt$assoNew
      power2 <- sparsReslt$power
      simMat2 <- .transToSim(x = dissMat2, simFunc = simFunc, 
                             simFuncPar = simFuncPar)
      adjaMat2 <- .transToAdja(x = simMat2, weighted = weighted)
    }
    else {
      dissEst2 <- dissScale2 <- assoMat2 <- assoEst2 <- dissMat2 <- power2 <- simMat2 <- adjaMat2 <- NULL
    }
    assoBoot1 <- assoBoot2 <- NULL
  }
  else {
    if (verbose %in% 2:3) {
      if (sparsMethod != "none") {
        message("\nSparsify associations via '", sparsMethod, 
                "' ... ", appendLF = FALSE)
      }
    }
    sparsReslt <- .sparsify(assoMat = assoMat1, countMat = counts1, 
                            sampleSize = sampleSize[1], measure = measure, measurePar = measurePar, 
                            assoType = assoType, sparsMethod = sparsMethod, thresh = thresh[1], 
                            alpha = alpha[1], adjust = adjust, lfdrThresh = lfdrThresh[1], 
                            trueNullMethod = trueNullMethod, nboot = nboot, assoBoot = assoBoot, 
                            softThreshType = softThreshType, softThreshPower = softThreshPower[1], 
                            softThreshCut = softThreshCut[1], cores = cores, 
                            logFile = logFile, kNeighbor = kNeighbor, knnMutual = knnMutual, 
                            verbose = verbose, seed = seed)
    if (verbose %in% 2:3 & sparsMethod != "none") 
      message("Done.")
    assoEst1 <- assoMat1
    assoMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    dissEst1 <- dissScale1 <- NULL
    dissMat1 <- .transToDiss(x = assoMat1, dissFunc = dissFunc, 
                             dissFuncPar = dissFuncPar)
    if (sparsMethod == "softThreshold") {
      simMat1 <- sparsReslt$simMat
      adjaMat1 <- assoMat1
      assoBoot1 <- NULL
    }
    else {
      simMat1 <- .transToSim(x = dissMat1, simFunc = simFunc, 
                             simFuncPar = simFuncPar)
      adjaMat1 <- .transToAdja(x = simMat1, weighted = weighted)
      if (sparsMethod == "bootstrap" && !is.null(assoBoot) && 
          !is.list(assoBoot) && assoBoot == TRUE) {
        assoBoot1 <- sparsReslt$assoBoot
      }
      else {
        assoBoot1 <- NULL
      }
    }
    if (twoNets) {
      if (verbose %in% 2:3) {
        if (sparsMethod != "none") {
          message("\nSparsify associations in group 2 ... ", 
                  appendLF = FALSE)
        }
      }
      sparsReslt <- .sparsify(assoMat = assoMat2, countMat = counts2, 
                              sampleSize = sampleSize[2], measure = measure, 
                              measurePar = measurePar, assoType = assoType, 
                              sparsMethod = sparsMethod, thresh = thresh[2], 
                              alpha = alpha[2], adjust = adjust, lfdrThresh = lfdrThresh[2], 
                              trueNullMethod = trueNullMethod, nboot = nboot, 
                              assoBoot = assoBoot, softThreshType = softThreshType, 
                              softThreshPower = softThreshPower[2], softThreshCut = softThreshCut[2], 
                              cores = cores, logFile = logFile, kNeighbor = kNeighbor, 
                              knnMutual = knnMutual, verbose = verbose, seed = seed)
      if (verbose %in% 2:3 & sparsMethod != "none") 
        message("Done.")
      assoEst2 <- assoMat2
      assoMat2 <- sparsReslt$assoNew
      power2 <- sparsReslt$power
      dissEst2 <- dissScale2 <- NULL
      dissMat2 <- .transToDiss(x = assoMat2, dissFunc = dissFunc, 
                               dissFuncPar = dissFuncPar)
      if (sparsMethod == "softThreshold") {
        simMat2 <- sparsReslt$simMat
        adjaMat2 <- assoMat2
        assoBoot2 <- NULL
      }
      else {
        simMat2 <- .transToSim(x = dissMat2, simFunc = simFunc, 
                               simFuncPar = simFuncPar)
        adjaMat2 <- .transToAdja(x = simMat2, weighted = weighted)
        if (sparsMethod == "bootstrap" && !is.null(assoBoot) && 
            !is.list(assoBoot) && assoBoot == TRUE) {
          assoBoot2 <- sparsReslt$assoBoot
        }
        else {
          assoBoot2 <- NULL
        }
      }
    }
    else {
      dissEst2 <- dissScale2 <- assoMat2 <- assoEst2 <- dissMat2 <- power2 <- simMat2 <- adjaMat2 <- assoBoot2 <- NULL
    }
  }
  g <- graph_from_adjacency_matrix(adjaMat1, weighted = TRUE, 
                                   mode = "undirected", diag = FALSE)
  if (is.null(E(g)$weight)) {
    isempty1 <- TRUE
    edgelist1 <- NULL
  }
  else {
    isempty1 <- FALSE
    edgelist1 <- data.frame(get.edgelist(g))
    colnames(edgelist1) <- c("v1", "v2")
    if (!is.null(assoMat1)) {
      edgelist1$asso <- sapply(1:nrow(edgelist1), function(i) {
        assoMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })
    }
    edgelist1$diss <- sapply(1:nrow(edgelist1), function(i) {
      dissMat1[edgelist1[i, 1], edgelist1[i, 2]]
    })
    if (all(adjaMat1 %in% c(0, 1))) {
      edgelist1$sim <- sapply(1:nrow(edgelist1), function(i) {
        simMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })
    }
    edgelist1$adja <- sapply(1:nrow(edgelist1), function(i) {
      adjaMat1[edgelist1[i, 1], edgelist1[i, 2]]
    })
  }
  if (twoNets) {
    g <- graph_from_adjacency_matrix(adjaMat2, weighted = TRUE, 
                                     mode = "undirected", diag = FALSE)
    if (is.null(E(g)$weight)) {
      isempty2 <- TRUE
      edgelist2 <- NULL
    }
    else {
      isempty2 <- FALSE
      edgelist2 <- data.frame(get.edgelist(g))
      colnames(edgelist2) <- c("v1", "v2")
      if (!is.null(assoMat2)) {
        edgelist2$asso <- sapply(1:nrow(edgelist2), function(i) {
          assoMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }
      edgelist2$diss <- sapply(1:nrow(edgelist2), function(i) {
        dissMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })
      if (all(adjaMat2 %in% c(0, 1))) {
        edgelist2$sim <- sapply(1:nrow(edgelist2), function(i) {
          simMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }
      edgelist2$adja <- sapply(1:nrow(edgelist2), function(i) {
        adjaMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })
    }
    if (isempty1 && verbose > 0) {
      message("\nNetwork 1 has no edges.")
    }
    if (isempty2 && verbose > 0) {
      message("Network 2 has no edges.")
    }
  }
  else {
    edgelist2 <- NULL
    if (isempty1 && verbose > 0) {
      message("\nNetwork has no edges.")
    }
  }
  output <- list()
  output$edgelist1 <- edgelist1
  output$edgelist2 <- edgelist2
  output$assoMat1 <- assoMat1
  output$assoMat2 <- assoMat2
  output$dissMat1 <- dissMat1
  output$dissMat2 <- dissMat2
  output$simMat1 <- simMat1
  output$simMat2 <- simMat2
  output$adjaMat1 <- adjaMat1
  output$adjaMat2 <- adjaMat2
  output$assoEst1 <- assoEst1
  output$assoEst2 <- assoEst2
  output$dissEst1 <- dissEst1
  output$dissEst2 <- dissEst2
  output$dissScale1 <- dissScale1
  output$dissScale2 <- dissScale2
  output$assoBoot1 <- assoBoot1
  output$assoBoot2 <- assoBoot2
  output$countMat1 <- countsOrig1
  output$countMat2 <- countsOrig2
  if (!is.null(countsJointOrig)) 
    output$countsJoint <- countsJointOrig
  output$normCounts1 <- counts1
  output$normCounts2 <- counts2
  output$groups <- groups
  output$matchDesign <- matchDesign
  output$sampleSize <- sampleSize
  output$softThreshPower <- list(power1 = power1, power2 = power2)
  output$assoType <- assoType
  output$twoNets <- twoNets
  output$parameters <- list(dataType = dataType, group = group, 
                            filtTax = filtTax, filtTaxPar = filtTaxPar, filtSamp = filtSamp, 
                            filtSampPar = filtSampPar, jointPrepro = jointPrepro, 
                            zeroMethod = zeroMethod, zeroPar = zeroPar, needfrac = needfrac, 
                            needint = needint, normMethod = normMethod, normPar = normPar, 
                            measure = measure, measurePar = measurePar, sparsMethod = sparsMethod, 
                            thresh = thresh, adjust = adjust, trueNullMethod = trueNullMethod, 
                            alpha = alpha, lfdrThresh = lfdrThresh, nboot = nboot, 
                            softThreshType = softThreshType, softThreshPower = softThreshPower, 
                            softThreshCut = softThreshCut, kNeighbor = kNeighbor, 
                            knnMutual = knnMutual, dissFunc = dissFunc, dissFuncPar = dissFuncPar, 
                            simFunc = simFunc, simFuncPar = simFuncPar, scaleDiss = scaleDiss, 
                            weighted = weighted, sampleSize = sampleSize)
  output$call = match.call()
  class(output) <- "microNet"
  return(output)
}
<bytecode: 0x000001d1d171de10>
  <environment: namespace:NetCoMi>