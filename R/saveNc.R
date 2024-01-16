saveNetCdf <- function(outNcFilePath, name, dimList, values) {
  outNc <- create.nc(outNcFilePath, format = "netcdf4", share = FALSE)
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
  var.def.nc(outNc, name, "NC_DOUBLE", names(dimList), deflate = 9)
  var.put.nc(outNc, name, values)
  close.nc(outNc)
}

