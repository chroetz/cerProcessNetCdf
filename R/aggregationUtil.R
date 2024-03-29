.info <- new.env()


clearInfo <- function() {
  rm(list=ls(.info), envir=.info)
}


assertLonLat <- function(increasingLonValues, increasingLatValues) {
  stopifnot(
    max(abs(.info$grid$increasingLatValues - increasingLatValues)) < .info$gridTol,
    max(abs(.info$grid$increasingLonValues - increasingLonValues)) < .info$gridTol)
}


getSingleBoundingBox <- function(boundingBoxes, regionName) {
  if (is.null(boundingBoxes)) {
    return(NULL)
  }
  boundingBoxIndex <- which(boundingBoxes$regionName == regionName)
  stopifnot(length(boundingBoxIndex) == 1)
  nms <- c("min_lon", "max_lon", "min_lat", "max_lat")
  bbInfo <- lapply(nms, \(nm) boundingBoxes[[nm]][boundingBoxIndex])
  names(bbInfo) <- nms
  return(bbInfo)
}


subsetBox <- function(valuesOnTotalGrid, bbInfo) {
  # TODO: data source format and bounding box format
  res <- valuesOnTotalGrid[
    bbInfo$min_lon:bbInfo$max_lon,
    bbInfo$min_lat:bbInfo$max_lat,
    drop=FALSE]
  return(res)
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
  nc <- openNc(filePath)
  maskList$nc <- nc
  maskList$lonValues <- var.get.nc(nc, "lon") |> as.vector()
  maskList$latValues <- var.get.nc(nc, "lat") |> as.vector()
  assertLonLat(maskList$lonValues, rev(maskList$latValues))
  maskList$lonLatIdx <-
    c(lon = ncGetDimensionIndex(nc, "lon"),
      lat = ncGetDimensionIndex(nc, "lat"))
  maskList$regionNames <- ncGetNonDimVariableNames(nc)
  maskList$gridFormat <- getNativeGridFormatFromNc(nc)
  cat("Grid format of mask file:", format(maskList$gridFormat), "\n")
  stopifnot(identical(.info$targetFormat, maskList$gridFormat))
  .info$maskList <- maskList
}


getMaskValues <- function(regionName, maskList, bbInfo = NULL) {
  # TODO: respect target format
  if (hasValue(bbInfo)) {
    # TODO: assumes lon, lat order and same format as bounding box
    values <- var.get.nc(
      maskList$nc,
      regionName,
      start = c(bbInfo$min_lon, bbInfo$min_lat),
      count = c(bbInfo$max_lon - bbInfo$min_lon + 1, bbInfo$max_lat - bbInfo$min_lat + 1),
      collapse = FALSE)
  } else {
    values <- var.get.nc(
      maskList$nc,
      regionName,
      collapse = FALSE)
  }

  if (any(is.na(values))) {
    message("WARNING: NAs in mask values in region ", regionName)
  }

  return(values)
}


readBoundingBoxes <- function(filePath) {
  nc <- openNc(filePath)
  on.exit(close.nc(nc))
  bb <- read.nc(nc)
  bbFormat <- getNativeGridFormatFromNc(nc, onlyLonLat=TRUE)

  # TODO: this is a bit hacky:
  # The region variable is the only character variable. Make sure it is named "regionName"
  selRegion <- sapply(bb, typeof) == "character"
  names(bb)[selRegion] <- "regionName"

  cat(
    "Grid format of bounding boxes:",
    format(bbFormat),
    "\n")
  stopifnot(
    bbFormat$nLon == .info$maskList$gridFormat$nLon,
    bbFormat$nLat == .info$maskList$gridFormat$nLat,
    bbFormat$lonIncreasing == .info$maskList$gridFormat$lonIncreasing,
    bbFormat$latIncreasing == .info$maskList$gridFormat$latIncreasing)
  .info$boundingBoxes <- bb
  .info$boundingBoxFormat <- bbFormat
}


readMaskSum <- function(filePath) {
  nc <- openNc(filePath)
  on.exit(close.nc(nc))

  gridFormat <- getNativeGridFormatFromNc(nc)
  cat(
    "Grid format of mask sum:",
    format(gridFormat),
    "\n")
  maskSum <- var.get.nc(nc, "maskSum")

  .info$maskSum <- ensureGridFormat(
    maskSum,
    .info$maskList$gridFormat,
    gridFormat)
}

