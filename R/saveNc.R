saveNetCdf <- function(outFilePath, dimList, valueList, deflate = 9) {
  outNc <- create.nc(outFilePath, format = "netcdf4", share = FALSE)
  for (i in seq_along(dimList)) {
    dimName <- names(dimList)[i]
    dimType <- inferNetCdfType(dimList[[i]])
    dim.def.nc(outNc, dimName, dimlength = length(dimList[[i]]))
    var.def.nc(outNc, dimName, dimType, dimName)
    var.put.nc(outNc, dimName, dimList[[i]])
    attris <- attributes(dimList[[i]])
    for (nm in names(attris)) {
      att.put.nc(outNc, dimName, nm, inferNetCdfType(attris[[nm]]), attris[[nm]])
    }
  }
  for (nm in names(valueList)) {
    var.def.nc(outNc, nm, "NC_DOUBLE", names(dimList), deflate = deflate)
    var.put.nc(outNc, nm, valueList[[nm]])
  }
  close.nc(outNc)
}


initNetCdf <- function(outFilePath, dimList, varNames, deflate = 9) {
  outNc <- create.nc(outFilePath, format = "netcdf4", share = FALSE)
  for (i in seq_along(dimList)) {
    dimName <- names(dimList)[i]
    dimType <- inferNetCdfType(dimList[[i]])
    dim.def.nc(outNc, dimName, dimlength = length(dimList[[i]]))
    var.def.nc(outNc, dimName, dimType, dimName)
    var.put.nc(outNc, dimName, dimList[[i]])
    attris <- attributes(dimList[[i]])
    for (nm in names(attris)) {
      att.put.nc(outNc, dimName, nm, inferNetCdfType(attris[[nm]]), attris[[nm]])
    }
  }
  for (nm in varNames) {
    var.def.nc(outNc, nm, "NC_DOUBLE", names(dimList), deflate = deflate)
  }
  close.nc(outNc)
}


initCopyNetCdf <- function(outFilePath, sourceFilePath, deflate = 9, ...) {

  overwriteDim <- list(...)

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
      if (varInq$name %in% names(overwriteDim)) {
        var.put.nc(outNc, varInq$name, overwriteDim[[varInq$name]])
      } else {
        var.put.nc(outNc, varInq$name, var.get.nc(nc, varInq$name))
      }
    }
    if (varInq$name %in% names(overwriteDim)) {
      attris <- attributes(overwriteDim[[varInq$name]])
      for (nm in names(attris)) {
        att.put.nc(
          outNc,
          varInq$name,
          nm,
          inferNetCdfType(attris[[nm]]),
          attris[[nm]])
      }
    } else {
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
  }

  close.nc(nc)

  cat("done after", (proc.time()-pt)[3], "s\n")

  return(outNc)
}


saveLonLatTimeToNetCdf <- function(nc, info, data) {
  lonIdx <- which(info$lon %in% dimnames(data)[["lon"]])
  latIdx <- which(info$lat %in% dimnames(data)[["lat"]])
  timeIdx <- which(info$timeInterpreted %in% dimnames(data)[["time"]])
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





inferNetCdfType <- function(x) {
  switch(
    typeof(x),
    "integer" = "NC_INT",
    "double" = "NC_DOUBLE",
    "character" = "NC_STRING",
    stop("Cannot process type: ", typeof(x)))
}
