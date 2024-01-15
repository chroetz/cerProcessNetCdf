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

  cat("Start\n")
  batchIndices <- seq_len(nBatches)
  if (hasValue(batchIndexFilter)) {
    batchIndices <- intersect(batchIndices, batchIndexFilter)
  }
  for (k in batchIndices) {

    # Don't overwrite
    if (file.exists(outFileNames[k])) next

    indicesOfBatch <- listOfIndicesOfBatches[[k]]

    cat("Start Batch ", k, "/", nBatches,
        "from", min(indicesOfBatch), "to", max(indicesOfBatch), "\n")

    pt <- proc.time()[3]

    maskArray <- sapply(
      indicesOfBatch,
      \(index) {
        ptInner <- proc.time()[3]
        cat(index, "... ")
        res <- exactextractr::coverage_fraction(globe$raster, sf[index, ])
        # In the following line, t() turns lat lon order to lon lat, which is the format assumed by writeMasksAsNetCdf().
        mat <- t(raster::as.matrix(res[[1]]))
        cat("done after", proc.time()[3] - ptInner, "s\n")
        return(mat)
      },
      simplify="array")

    cat("\n")
    cat("Finished after", proc.time()[3] - pt, "s\n")

    dimnames(maskArray) <- list(
      lon = character(0),
      lat = character(0),
      varName = sf[[idColumnName]][indicesOfBatch])

    writeMasksAsNetCdf(
      outFileNames[k],
      maskArray,
      dimLon = globe$dimLon, dimLat = globe$dimLat)
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

  maskArray <- sapply(
    shapeFilePaths,
    simplify = "array",
    rasterizeShape,
    raster = globe$raster
  )

  stopifnot(dim(maskArray)[1:2] == c(nLon, nLat))

  dimnames(maskArray) <- list(
    lon = character(0),
    lat = character(0),
    region = regionNames)

  writeMasksAsNetCdf(
    outFilePath,
    maskArray,
    dimLon = globe$dimLon,
    dimLat = globe$dimLat)
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


writeMasksAsNetCdf <- function(outFilePath, maskArray, dimLon, dimLat) {

  pt <- proc.time()
  cat("Write mask data of size", dim(maskArray), "to", outFilePath, "... ")

  stopifnot(
    length(dimLon) == dim(maskArray)[1],
    length(dimLat) == dim(maskArray)[2])
  nLon <- length(dimLon)
  nLat <- length(dimLat)
  varNames <- dimnames(maskArray)[[3]]
  rnc <- create.nc(outFilePath, format = "netcdf4")
  dim.def.nc(rnc, "lon", dimlength = nLon)
  var.def.nc(rnc, "lon", "NC_DOUBLE", "lon")
  var.put.nc(rnc, "lon", dimLon)
  att.put.nc(rnc, "lon", "units", "NC_CHAR", "degree east")
  dim.def.nc(rnc, "lat", dimlength = nLat)
  var.def.nc(rnc, "lat", "NC_DOUBLE", "lat")
  var.put.nc(rnc, "lat", dimLat)
  att.put.nc(rnc, "lat", "units", "NC_CHAR", "degree north")
  for (varName in varNames) {
    var.def.nc(rnc, varName, "NC_DOUBLE", c("lon", "lat"), deflate = 9)
    att.put.nc(rnc, varName, "units", "NC_CHAR", "1")
    var.put.nc(rnc, varName, maskArray[,,varName])
  }
  close.nc(rnc)

  cat("done after", (proc.time() - pt)[3], "s\n")
}
