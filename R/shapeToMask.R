#' @export
runShapeToMaskOneFileForAllRegions <- function(
    shapeFilePath,
    nLon, nLat,
    outFilePath,
    metaOutFilePath = NULL,
    idColumnName = NULL,
    nBatches = 1,
    batchIndex = 1
) {

  cat("read geodata file...")
  pt <- proc.time()[3]
  sf <- sf::read_sf(shapeFilePath)
  cat(" done after", proc.time()[3] - pt, "s\n")

  isStandardClass <- sapply(sf, \(x) class(x)[1] %in% c("character", "numeric", "logical"))

  if (!hasValueString(idColumnName)) idColumnName <- guessIdColumnName(sf)

  if (hasValueString(metaOutFilePath)) {
    meta <-
      unclass(sf)[isStandardClass] |>
      as_tibble()
    readr::write_csv(meta, metaOutFilePath)
  }

  globe <- getGlobalRaster(nLon, nLat)

  # Prepare Batches.
  n <- nrow(sf)
  cat("Split regions into ", nBatches, "batches.\n")
  batches <- setupBatches(seq_len(n), nBatches)
  batch <- batches[[batchIndex]]
  if (length(batch) == 0) {
    cat("Batch", batchIndex, "is empty. Nothing to do.\n")
    return(invisible())
  }
  cat("Process batch", batchIndex, "with", length(batch), "regions\n")

  if (length(batch) == 0 | any(is.na(batch))) {
    cat("Empty or invalid batch. Skipping.\n")
    return(invisible())
  }

  cat("Start batch from", min(batch), "to", max(batch), "\n")

  if (file.exists(outFilePath)) {
    cat("Skip batch because", outFilePath, "already exists.\n")
    return(invisible())
  }

  pt <- proc.time()[3]

  outNc <- initLonLatNetCdf(outFilePath, globe)

  for (index in batch) {
    ptInner <- proc.time()[3]
    cat(index, "... ")
    res <- exactextractr::coverage_fraction(globe$raster, sf[index, ])
    mat <- t(raster::as.matrix(res[[1]]))
    varName <- sf[[idColumnName]][index]
    writeLonLatVariable(outNc, varName, mat)
    cat("done after", proc.time()[3] - ptInner, "s\n")
  }

  close.nc(outNc)
  cat("Finished after", proc.time()[3] - pt, "s\n")
}


#' @export
runShapeToMaskOneFilePerRegion <- function(
    shapeFilePaths,
    nLon, nLat,
    outFilePath
) {

  regionNames <- names(shapeFilePaths)
  if (!hasValue(regionNames)) {
    regionNames <- removeFileNameEnding(basename(shapeFilePaths))
  }
  if (length(unique(regionNames)) != length(regionNames)) {
    stop("Unable to infer unqiue names for regions.")
  }

  globe <- getGlobalRaster(nLon, nLat)

  outNc <- initLonLatNetCdf(outFilePath, globe)

  for (k in seq_along(shapeFilePaths)) {
    mat <- rasterizeShape(shapeFilePaths[k], globe$raster)
    stopifnot(dim(mat) == c(nLon, nLat))
    writeLonLatVariable(outNc, regionNames[k], mat)
  }

  close.nc(outNc)
}


rasterizeShape <- function(shapeFilePath, raster) {
  cat("process", shapeFilePath, "\n")
  regionShape <- sf::read_sf(shapeFilePath)
  mask <- exactextractr::coverage_fraction(raster, regionShape)
  stopifnot(length(mask) == 1)
  res <- t(raster::as.matrix(mask[[1]]))
  return(res)
}


guessIdColumnName <- function(sf) {
  isString <- sapply(sf, is.character)
  uniqueVals <- sapply(sf[isString], \(x) length(unique(x)))
  return(names(sf)[isString][which.max(uniqueVals)])
}


getGlobalRaster <- function(nLon, nLat) {

  rangeLon <- c(-180, 180)
  rangeLat <- c(-90, 90)

  rasterGlobal <- raster::raster(
    nrows = nLat,
    ncols = nLon,
    xmn = rangeLon[1], xmx = rangeLon[2],
    ymn = rangeLat[1], ymx = rangeLat[2])

  return(list(
    raster = rasterGlobal,
    dimLon = raster::xFromCol(rasterGlobal),
    dimLat = raster::yFromRow(rasterGlobal),
    coordinates = raster::coordinates(rasterGlobal)))
}


initLonLatNetCdf <- function(outFilePath, raster) {
  cat("Create output file", outFilePath, "\n")
  outNc <- create.nc(outFilePath, format = "netcdf4")
  dim.def.nc(outNc, "lon", dimlength = raster$dimLon |> length())
  var.def.nc(outNc, "lon", "NC_DOUBLE", "lon")
  var.put.nc(outNc, "lon", raster$dimLon)
  att.put.nc(outNc, "lon", "units", "NC_CHAR", "degree east")
  dim.def.nc(outNc, "lat", dimlength = raster$dimLat |> length())
  var.def.nc(outNc, "lat", "NC_DOUBLE", "lat")
  var.put.nc(outNc, "lat", raster$dimLat)
  att.put.nc(outNc, "lat", "units", "NC_CHAR", "degree north")
  return(outNc)
}


writeLonLatVariable <- function(outNc, varName, data, units = "1") {
  var.def.nc(outNc, varName, "NC_DOUBLE", c("lon", "lat"), deflate = 9)
  att.put.nc(outNc, varName, "units", "NC_CHAR", units)
  var.put.nc(outNc, varName, data)
}
