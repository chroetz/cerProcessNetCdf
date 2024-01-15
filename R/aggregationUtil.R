.info <- new.env()


processRegionYear <- function(regionName, year, invarNames, aggregationDistri, batchSize, outNc) {
  stopifnot(batchSize >= 1)
  n <- length(invarNames)
  nBatches <- ceiling(n/batchSize)
  batchIdxs <- lapply(
    seq_len(nBatches),
    \(k) {
      idxs <- (1+(k-1)*batchSize):(k*batchSize)
      idxs[idxs <= n]
    }
  )
  for (batchNr in seq_len(nBatches)) {
    idxs <- batchIdxs[[batchNr]]
    cat("\t\tVariable indices from", min(idxs), "to", max(idxs), "\n")
    pt <- proc.time()
    invarValues <- getInvarValues(year, min(idxs), length(idxs), regionName) # this takes time
    cat("\t\t\tgetInvarValues() took", (proc.time()-pt)[3], "s\n")
    pt <- proc.time()
    for (statisticName in .info$statisticNames) {
      cat("\t\t\tStatistic:", statisticName, "...")
      x <- calculateStatisticOnGrid(statisticName, invarValues)
      results <- integrateDistribution(aggregationDistri, x)
      saveResult(results, year, regionName, statisticName, invarNames[idxs], outNc=outNc)
      cat(" Done.\n")
    }
    cat("\t\t\tcalculating and saving took", (proc.time()-pt)[3], "s\n")
  }
}



getFullyFilledRegionNames <- function(info, year, invarNames, outNc) {
  regionNames <- var.get.nc(outNc, "region")
  variableNames <- ncGetNonDimVariableNames(outNc)

  if (!is.character(invarNames)) {
    invarNames <- paste0(info$invarValueVariableName, "_", invarNames)
  }

  if (any(!invarNames %in% variableNames)) {
    return(NULL)
  }
  var1 <- tryCatch(var.inq.nc(outNc, 1), error = \(cond) FALSE)
  if (isFALSE(var1)) {
    outNcFilePath <- getOutNcFilePath(year)
    stop("The file ", outNcFilePath, " is corrupt! Probably need to delete it and run calculations again.")
  }
  hasNa <- sapply(
    invarNames,
    \(invarName) rowSums(is.na(var.get.nc(outNc, invarName, collapse=FALSE))) > 0)
  isRegionFilled <- rowSums(hasNa) == 0
  return(regionNames[isRegionFilled])
}


assertLonLat <- function(increasingLonValues, increasingLatValues) {
  stopifnot(
    max(abs(.info$grid$increasingLatValues - increasingLatValues)) < .info$gridTol,
    max(abs(.info$grid$increasingLonValues - increasingLonValues)) < .info$gridTol)
}

calculateStatisticOnGrid <- function(statisticName, invarValues) {
  .info$statisticFunctions[[statisticName]](invarValues)
}



calculateProductDistribution <- function(factor1, factor2, naToZero = TRUE) {
  stopifnot(identical(dim(factor2), dim(factor1)))
  product <- factor2 * factor1
  distibution <- normalizeDistribution(product, naToZero = naToZero)
  return(distibution)
}


normalizeDistribution <- function(unnormalizedDistribution, naToZero = TRUE) {
  distibution <- unnormalizedDistribution / sum(unnormalizedDistribution, na.rm = TRUE)
  if (naToZero) {
    distibution[is.na(distibution)] <- 0
  }
  stopifnot(abs(sum(distibution) - 1) < sqrt(.Machine$double.eps))
  return(distibution)
}


getSingleBoundingBox <- function(boundingBoxes, regionName) {
  boundingBoxIndex <- which(boundingBoxes$GID_1 == regionName)
  stopifnot(length(boundingBoxIndex) == 1)
  nms <- c("min_lon", "max_lon", "min_lat", "max_lat")
  bbInfo <- lapply(nms, \(nm) boundingBoxes[[nm]][boundingBoxIndex])
  names(bbInfo) <- nms
  return(bbInfo)
}


subsetRegion <- function(valuesOnTotalGrid, regionName, boundingBoxes, grid) {
  bbInfo <- getSingleBoundingBox(boundingBoxes, regionName)
  len <- grid$latValues |> length()
  # TODO: check whether one needs to reverse the lat values
  reversedLat <- reverseIndex(bbInfo$min_lat, bbInfo$max_lat, len)
  res <- valuesOnTotalGrid[
    bbInfo$min_lon:bbInfo$max_lon,
    reversedLat$from:reversedLat$to,
    drop=FALSE]
  return(res)
}


getInvarFilePath <- function(year) {
  fileInfo <- .info$invarFileMeta |> filter(.data$year == .env$year)
  stopifnot(nrow(fileInfo) == 1)
  return(fileInfo$path)
}


getInvarNames <- function(year) {
  filePath <- getInvarFilePath(year)
  nc <- open.nc(filePath)
  invarNames <- var.get.nc(nc, .info$invarDimensionName)
  close.nc(nc)
  return(invarNames)
}


reverseArrayDim <- function(x, dims) {
  if (is.character(dims)) {
    dims <- vapply(dims, \(nm) which(nm == names(dimnames(x)))[[1]], integer(1))
  }
  dims <- as.integer(dims)
  n <- length(dim(x))
  stopifnot(all(dims > 0 & dims <= n))
  dimIndices <- rep(list(rlang::missing_arg()), n)
  for (i in dims) {
    dimLen <- dim(x)[i]
    if (dimLen == 0) next
    dimIndices[[i]] <- dimLen:1
  }
  z <- do.call(`[`, c(list(x), dimIndices))
  oldClass(z) <- oldClass(x)
  return(z)
}


reverseIndex <- function(from, to, len) {
  list(
    from = len - to + 1,
    to = len - from + 1,
    count = abs(from - to) + 1)
}


openAndCheckMaskFile <- function(filePath) {
  maskList <- list()
  nc <- open.nc(filePath)
  maskList$nc <- nc
  maskList$lonValues <- var.get.nc(nc, "lon")
  maskList$latValues <- var.get.nc(nc, "lat")
  assertLonLat(maskList$lonValues, rev(maskList$latValues))
  maskList$lonLatIdx <-
    c(lon = ncGetDimensionIndex(nc, "lon"),
      lat = ncGetDimensionIndex(nc, "lat"))
  maskList$regionNames <- ncGetNonDimVariableNames(nc)
  .info$maskList <- maskList
}


getMaskValues <- function(regionName, maskList, boundingBoxes = NULL) {
  if (hasValue(boundingBoxes)) {
    bbInfo <- getSingleBoundingBox(boundingBoxes, regionName)
    values <- var.get.nc(
      maskList$nc,
      regionName,
      start = c(bbInfo$min_lon, bbInfo$min_lat),
      count = c(bbInfo$max_lon - bbInfo$min_lon + 1, bbInfo$max_lat - bbInfo$min_lat + 1),
      collapse = FALSE)
  } else {
    values <- var.get.nc(maskList$nc, regionName)
  }

  values <- reverseArrayDim(values, maskList$lonLatIdx["lat"])

  if (any(is.na(values))) {
    message("WARNING: NAs in mask values in region ", regionName)
  }

  return(values)
}


checkInvar <- function(year, invarNames) {
  filePath <- getInvarFilePath(year)
  invar <- list()
  nc <- open.nc(filePath)
  invar$lonValues <- var.get.nc(nc, "lon")
  invar$latValues <- var.get.nc(nc, "lat")
  invar$dimensionValues <- var.get.nc(nc, .info$invarDimensionName)
  close.nc(nc)
  assertLonLat(invar$lonValues, invar$latValues)
  stopifnot(all(invarNames %in% invar$dimensionValues))
  return(invisible())
}


getInvarValues <- function(year, fromIdx, count, regionName) {
  bbInfo <- .info$idxBoundingBoxes |> filter(.data$GID_1 == regionName) |> as.list()
  filePath <- getInvarFilePath(year)
  # lat is reversed in the mask, i.e., also in the bounding box
  nc <- open.nc(filePath)
  latCount <- bbInfo$max_lat - bbInfo$min_lat + 1
  lonCount <- bbInfo$max_lon - bbInfo$min_lon + 1
  len <- .info$grid$latValues |> length()
  reversedLat <- reverseIndex(bbInfo$min_lat, bbInfo$max_lat, len)
  invarValues <- var.get.nc(
    nc,
    .info$invarValueVariableName,
    start = c(bbInfo$min_lon, reversedLat$from, fromIdx)[.info$lonLatVarToDimOrder],
    count = c(lonCount, latCount, count)[.info$lonLatVarToDimOrder],
    collapse = FALSE)
  close.nc(nc)
  if (any(is.na(invarValues))) {
    message("WARNING: NAs in invar values in year ", year, ", `fromIdx` ", fromIdx, ", `count` ", count)
  }
  if (!all(.info$lonLatVarToDimOrder == 1:3)) {
    invarValues <- aperm(invarValues, order(.info$lonLatVarToDimOrder))
  }
  stopifnot(length(dim(invarValues)) == 3)
  stopifnot(dim(invarValues)[3] == count)
  return(invarValues)
}


integrateDistribution <- function(distri, values) {
  stopifnot(identical(dim(distri), dim(values)[1:2]))
  stopifnot(length(dim(values)) == 3)
  integral <- apply(values, 3, \(x) sum(distri * x))
  return(integral)
}


getOutNcFilePath <- function(year) {
  file.path(
    .info$outDir,
    sprintf(.info$outNcFilePattern, year))
}


initOutNc <- function(years, regionNames, statisticNames) {
  for (year in years) {
    outNcFilePath <- getOutNcFilePath(year)
    if (file.exists(outNcFilePath)) {
      cat(outNcFilePath, " already exits. Do not recreate.\n")
      next
    }
    outNc <- create.nc(outNcFilePath, format = "netcdf4", share = FALSE)
    dim.def.nc(outNc, "region", dimlength = length(regionNames))
    var.def.nc(outNc, "region", "NC_STRING", "region")
    var.put.nc(outNc, "region", regionNames)
    dim.def.nc(outNc, "statistic", dimlength = length(statisticNames))
    var.def.nc(outNc, "statistic", "NC_STRING", "statistic")
    var.put.nc(outNc, "statistic", statisticNames)
    close.nc(outNc)
  }
}




saveResult <- function(results, year, regionName, statisticName, variableNames, outNc) {
  stopifnot(length(results) == length(variableNames))
  if (any(is.na(results))) {
    message(
      "Got NA results in year ", year,
      ", region ", regionName,
      ", statistic ", statisticName,
      ", variables ", paste(variableNames[is.na(results)], collapse=", "))
  }
  regionNames <- var.get.nc(outNc, "region")
  regionIdx <- which(regionName == regionNames)
  stopifnot(length(regionIdx) == 1)
  statisticNames <- var.get.nc(outNc, "statistic")
  statisticIdx <- which(statisticName == statisticNames)
  stopifnot(length(statisticIdx) == 1)
  for (i in seq_along(variableNames)) {
    variableName <- variableNames[[i]]
    if (!is.character(variableName)) {
      variableName <- paste0(.info$invarValueVariableName, "_", variableName)
    }
    result <- results[[i]]
    if (!ncHasVariable(outNc, variableName)) {
      var.def.nc(outNc, variableName, "NC_DOUBLE", c("region", "statistic"), deflate = 9)
    }
    var.put.nc(
      outNc,
      variableName,
      result,
      start = c(regionIdx, statisticIdx),
      count = c(1, 1))
  }
}



readBoundingBoxes <- function(filePath) {
  nc <- open.nc(filePath)
  on.exit(close.nc(nc))
  .info$boundingBoxes <- read.nc(nc)
  # TODO: check GridFormat
}


readMaskSum <- function(filePath) {
  nc <- open.nc(filePath)
  on.exit(close.nc(nc))
  .info$maskSum <- read.nc(nc)
  # TODO: check GridFormat
}

