#' @export
concatAfterRoll <- function(
  dirPath,
  pattern,
  outDirPath,
  suffix = "",
  referenceDirPath,
  referencePattern,
  chunks = 10,
  deflate = 9
) {

  filePaths <- list.files(
    dirPath,
    pattern = pattern,
    full.names = TRUE)

  referenceFilePaths <- list.files(
    referenceDirPath,
    pattern = referencePattern,
    full.names = TRUE)

  meta <- getNetCdfDimensionMeta(filePaths)
  referenceMeta <- getNetCdfDimensionMeta(referenceFilePaths)

  for (rm in referenceMeta) {

    cat("Processing reference file", rm$filePath, "\n")

    outFilePath <- file.path(
      outDirPath,
      paste0(cerUtility::removeFileNameEnding(basename(rm$filePath)), suffix, ".nc"))

    outTimeDim <- rm$timeInterpreted |> as.integer()
    attributes(outTimeDim) <- list(units = "days since 1970-01-01", calendar = "gregorian")
    outNc <- initCopyNetCdf(outFilePath, rm$filePath, deflate = deflate, time = outTimeDim)
    # TODO: use initNetCdf() instead

    timeIdxAll <- seq_along(rm$timeInterpreted)
    timeIdxList <- cerUtility::setupBatches(timeIdxAll, chunks)

    for (timeIdxs in timeIdxList) {
      timeValues <- rm$timeInterpreted[timeIdxs]
      cat("\tGetting values for time value form ", min(timeValues), "to", max(timeValues), "... ")
      pt <- proc.time()
      v <- getAllForTimes(meta, timeValues)
      cat("done after", (proc.time()-pt)[3], "s\n")

      cat("\tSaving data to file ... ")
      pt <- proc.time()
      saveLonLatTimeToNetCdf(outNc, rm, v)
      cat("done after", (proc.time()-pt)[3], "s\n")
    }
    close.nc(outNc)
  }

  return(invisible(NULL))
}

getNetCdfDimensionMeta <- function(filePaths) {
  meta <- lapply(filePaths, \(filePath) {
    nc <- open.nc(filePath)
    res <- list(
      filePath = filePath,
      varName = ncGetNonDimVariableNames(nc),
      lon = var.get.nc(nc, "lon"),
      lat = var.get.nc(nc, "lat"),
      time = var.get.nc(nc, "time"),
      timeInterpreted = getNcTimes(nc, "time"))
    close.nc(nc)
    return(res)
  })
  return(meta)
}

getAllForTimes <- function(meta, timeValues) {
  dataList <- lapply(meta, \(m) {
    timeIdx <- which(m$timeInterpreted %in% timeValues)
    stopifnot(all(abs(diff(timeIdx)) == 1))
    if (length(timeIdx) != length(timeValues)) {
      stop("Some time values from the reference files are not available in the data files.")
    }
    start <- c(1, 1, min(timeIdx))
    count <- c(NA, NA, length(timeIdx))
    nc <- open.nc(m$filePath)
    data <- var.get.nc(
      nc,
      m$varName,
      start = start,
      count = count,
      collapse = FALSE)
    close.nc(nc)
    dimnames(data) <- list(lon = m$lon, lat = m$lat, time = timeValues)
    return(data)
  })
  values <- do.call(abind::abind, c(dataList, list(along=2, use.dnns=TRUE)))
  return(values)
}

