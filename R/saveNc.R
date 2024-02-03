saveNetCdf <- function(outFilePath, dimList, valueList) {
  outNc <- create.nc(outFilePath, format = "netcdf4", share = FALSE)
  for (i in seq_along(dimList)) {
    dimName <- names(dimList)[i]
    dimType <- switch(
      typeof(dimList[[i]]),
      "integer" = "NC_INT",
      "double" = "NC_DOUBLE",
      "character" = "NC_STRING",
      stop("Cannot process dimension type: ", typeof(dimList[[i]])))
    dim.def.nc(outNc, dimName, dimlength = length(dimList[[i]]))
    var.def.nc(outNc, dimName, dimType, dimName)
    var.put.nc(outNc, dimName, dimList[[i]])
  }
  for (nm in names(valueList)) {
    var.def.nc(outNc, nm, "NC_DOUBLE", names(dimList), deflate = 9)
    var.put.nc(outNc, nm, valueList[[nm]])
  }
  close.nc(outNc)
}

