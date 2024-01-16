#' GridFormat class
#' @description A class to store the format of a longitude - latitude grid of the globe.
#' @param nLon,nLat Number of longitude and latitude values.
#' @param lonIncreasing,latIncreasing Whether the longitude/latitude values are in increasing order.
#' @param lonFirst Whether the longitude is the first dimension.
#' @examples
#' createGridFormat(20, 10, TRUE, TRUE, TRUE)
#' @export
createGridFormat <- function(
    nLon,
    nLat,
    lonIncreasing,
    latIncreasing,
    lonFirst
) {
  gridFormat <- lst(
    nLon,
    nLat,
    lonIncreasing,
    latIncreasing,
    lonFirst)
  class(gridFormat) <- "GridFormat"
  return(gridFormat)
}

getLonValues <- function(n) {
  v <- getEquidistantCenters(n, -180, 180)
  return(v)
}

getLatValues <- function(n) {
  v <- getEquidistantCenters(n, -90, 90)
  return(v)
}

getEquidistantCenters <- function(n, from, to) {
  v <- seq(from, to, length.out = n+1)[-1] - (to-from)/(2*n)
  return(v)
}

#' @export
format.GridFormat <- function(x, ...) {
  paste0(
    "GridFormat:\n",
    "\tlon step size ", 360 / x$nLon, " degrees (", 360 / (x$nLon / 60)," arcmin)\n",
    "\tlat step size ", 180 / x$nLat, " degrees (", 180 / (x$nLat / 60)," arcmin)\n",
    "\tlon from -180 to 180 (", x$nLon, " values)\n",
    "\tlat from -90 to 90 (", x$nLat, " values)\n",
    "\torder: ", if (is.na(x$lonFirst)) "undefined" else if (x$lonFirst) "lon x lat" else "lat x lon", "\n",
    "\tlon in ", if (x$lonIncreasing) "increasing" else "decreasing", " order\n",
    "\tlat in ", if (x$latIncreasing) "increasing" else "decreasing", " order")
}


#' @export
print.GridFormat <- function(x, ...) {
  cat(format(x), "\n")
}


initializeGrid <- function(targetFormat) {

  cat("Initialize grid with target format:", format(targetFormat), "\n")

  .info$gridTol <- min(360 / targetFormat$nLon, 180 / targetFormat$nLat) / 10
  .info$grid <- c(
    targetFormat,
    list(
      increasingLonValues = getLonValues(targetFormat$nLon),
      increasingLatValues = getLatValues(targetFormat$nLat),
      orderedNames = if (targetFormat$lonFirst) c("lon", "lat") else c("lat", "lon")
    )
  )
  if (targetFormat$lonIncreasing) {
    .info$grid$lonValues <- .info$grid$increasingLonValues
  } else {
    .info$grid$lonValues <- .info$grid$increasingLonValues |> rev()
  }
  if (targetFormat$latIncreasing) {
    .info$grid$latValues <- .info$grid$increasingLatValues
  } else {
    .info$grid$latValues <- .info$grid$increasingLatValues |> rev()
  }
}


#' @export
getNativeGridFormatFromFile <- function(filePath, variableName = NULL, onlyLonLat = FALSE) {
  nc <- open.nc(filePath)
  on.exit(close.nc(nc))
  return(getNativeGridFormatFromNc(nc, variableName, onlyLonLat))
}


#' @export
getNativeGridFormatFromNc <- function(nc, variableName = NULL, onlyLonLat = FALSE) {

  allVarNames <- ncGetVariableNames(nc)
  stopifnot(c("lon", "lat") %in% allVarNames)
  lonValues <- var.get.nc(nc, "lon")
  latValues <- var.get.nc(nc, "lat")
  if (onlyLonLat) {
    lonFirst <- NA
  } else {
    if (!hasValue(variableName)) {
      variableName <- ncGetNonDimVariableNames(nc) |> first()
    }
    varInfo <- var.inq.nc(nc, variableName)
    stopifnot(varInfo$ndims >= 2)
    varDimNames <- vapply(
      varInfo$dimids,
      \(dimId) dim.inq.nc(nc, dimId)$name,
      character(1))
    stopifnot(c("lon", "lat") %in% varDimNames)

    lonFirst <- which(varDimNames == "lon") < which(varDimNames == "lat")
  }

  if (all(diff(lonValues) > 0)) {
    lonIncreasing <- TRUE
  } else if (all(diff(lonValues) < 0)) {
    lonIncreasing <- FALSE
  } else stop("Longitude values are neither increasing nor decreasing.")
  if (all(diff(latValues) > 0)) {
    latIncreasing <- TRUE
  } else if (all(diff(latValues) < 0)) {
    latIncreasing <- FALSE
  } else stop("Latitude values are neither increasing nor decreasing.")

  # check grid is uniform
  lonSteps <- abs(diff(lonValues))
  meanLonStep <- mean(lonSteps)
  stopifnot(max(abs(lonSteps - meanLonStep)) / meanLonStep < 0.1)
  stopifnot(min(lonValues) - meanLonStep < -180)
  stopifnot(max(lonValues) + meanLonStep > 180)
  latSteps <- abs(diff(latValues))
  meanLatStep <- mean(latSteps)
  stopifnot(max(abs(latSteps - meanLatStep)) / meanLatStep < 0.1)
  stopifnot(min(latValues) - meanLatStep < -90)
  stopifnot(max(latValues) + meanLatStep > 90)

  return(
    createGridFormat(
      nLon = length(lonValues),
      nLat = length(latValues),
      lonIncreasing = lonIncreasing,
      latIncreasing = latIncreasing,
      lonFirst = lonFirst))
}
