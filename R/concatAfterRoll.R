concatAfterRoll <- function(
  dirPath,
  pattern,
  outFilePath,
  referenceDirPath,
  referencePattern,
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

  for (i in seq_along(referenceMeta)) {
    rm <- referenceMeta[[i]]
    thisOutFilePath <- paste0(cerUtility::removeFileNameEnding(outFilePath), "_", i, ".nc")
    outNc <- initCopyNetCdf(thisOutFilePath, rm$filePath, deflate = deflate)
    timeValuesList <- list(rm$time[1:floor(length(rm$time)/2)], rm$time[(floor(length(rm$time)/2)+1):length(rm$time)])
    for (timeValues in timeValuesList) {
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
      time = var.get.nc(nc, "time"))
    close.nc(nc)
    return(res)
  })
  return(meta)
}

getAllForTimes <- function(meta, timeValues) {
  dataList <- lapply(meta, \(m) {
    timeIdx <- which(m$time %in% timeValues)
    stopifnot(all(abs(diff(timeIdx)) == 1))
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


initCopyNetCdf <- function(outFilePath, sourceFilePath, deflate = 9) {

  makeDirsIfNecessary(outFilePath)
  pt <- proc.time()
  cat("Creating output file", outFilePath, "... ")

  outNc <- create.nc(outFilePath, format = "netcdf4")
  nc <- open.nc(sourceFilePath)

  fileInq <- file.inq.nc(nc)
  dimNames <- ncGetDimensionNames(nc)

  for (idim in seq_len(fileInq$ndims)) {
    dimInq <- dim.inq.nc(nc, idim - 1)
    dim.def.nc(outNc, dimInq$name, dimInq$length, dimInq$unlim)
  }

  for (ivar in seq_len(fileInq$nvars)) {
    varInq <- var.inq.nc(nc, ivar - 1)
    defl <- if (varInq$name %in% dimNames) NA else deflate
    var.def.nc(outNc, varInq$name, varInq$type, varInq$dimids, deflate = defl)
    if (varInq$name %in% dimNames) {
      var.put.nc(outNc, varInq$name, var.get.nc(nc, varInq$name))
    }
    for (iatt in seq_len(varInq$natts)) {
      attInq <- att.inq.nc(nc, varInq$name, iatt - 1)
      att.put.nc(
        outNc,
        varInq$name,
        attInq$name,
        attInq$type,
        att.get.nc(nc, varInq$name, attInq$name))
    }
  }

  close.nc(nc)

  cat("done after", (proc.time()-pt)[3], "s\n")

  return(outNc)
}


saveLonLatTimeToNetCdf <- function(nc, info, data) {
  lonIdx <- which(info$lon == dimnames(data)[["lon"]])
  latIdx <- which(info$lat == dimnames(data)[["lat"]])
  timeIdx <- which(info$time == dimnames(data)[["time"]])
  start <- c(min(lonIdx), min(latIdx), min(timeIdx))
  count <- c(length(lonIdx), length(latIdx), length(timeIdx))
  var.put.nc(
    nc,
    info$varName,
    data,
    start = start,
    count = count)
  return(invisible(NULL))
}
