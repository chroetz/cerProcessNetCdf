#' @export
runShapeToMaskOneFileForAllRegions <- function(
    shapeFilePath,
    nLon, nLat,
    outFilePrefix,
    metaOutFilePath = NULL,
    idColumnName = NULL,
    nBatches = 10,
    batchIndexFilter = NULL
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
  allIndices <- seq_len(n)
  batchSize <- ceiling(n/nBatches)
  listOfIndicesOfBatches <- lapply(
    seq_len(nBatches),
    \(i) allIndices[((i-1)*batchSize+1):min(n, i*batchSize)])
  nDigits <- floor(log10(n)+1)
  outFileNames <- sapply(
    listOfIndicesOfBatches,
    \(indicesOfBatch) sprintf(
      paste0("%s_%0", nDigits,"d_to_%0", nDigits,"d.nc"),
      outFilePrefix, min(indicesOfBatch), max(indicesOfBatch)))
  batchIndices <- seq_len(nBatches)
  if (hasValue(batchIndexFilter)) {
    batchIndices <- intersect(batchIndices, batchIndexFilter)
  }
  if (length(batchIndices) == 0) {
    cat("No batches to process.\n")
    cat("nBatches:", nBatches, "\n")
    cat("batchSize:", batchSize, "\n")
    if (hasValue(batchIndexFilter)) {
      cat("batchIndexFilter:", paste(batchIndexFilter, collapse=","), "\n")
    }
    return(invisible())
  }

  cat("Process batches with index", paste(batchIndices, collapse=","), "\n")
  for (k in batchIndices) {

    indicesOfBatch <- listOfIndicesOfBatches[[k]]

    cat("Start Batch ", k, "/", nBatches,
        "from", min(indicesOfBatch), "to", max(indicesOfBatch), "\n")

    if (file.exists(outFileNames[k])) {
      cat("Skip batch", k, "because", outFileNames[k], "already exists.\n")
      next
    }

    pt <- proc.time()[3]

    outNc <- initLonLatNetCdf(outFileNames[k], globe)

    for (index in indicesOfBatch) {
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

  initLonLatNetCdf(outFilePath, globe)

  for (k in seq_along(shapeFilePaths)) {
    mat <- rasterizeShape(shapeFilePaths[k], globe$raster)
    stopifnot(dim(mat) == c(nLon, nLat))
    writeLonLatVariable(outFilePath, regionNames[k], mat)
  }
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
}


writeLonLatVariable <- function(outNc, varName, data, units = "1") {
  var.def.nc(outNc, varName, "NC_DOUBLE", c("lon", "lat"), deflate = 9)
  att.put.nc(outNc, varName, "units", "NC_CHAR", units)
  var.put.nc(outNc, varName, data)
}
