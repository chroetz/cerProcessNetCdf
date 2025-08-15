#' @export
ncHasVariable <- function(openNc, variableName) {
  ncInfo <- file.inq.nc(openNc)
  varInfo <- lapply(seq_len(ncInfo$nvars)-1, var.inq.nc, ncfile = openNc)
  varNames <- sapply(varInfo, \(vi) vi$name)
  variableName %in% varNames
}

#' @export
ncGetNonDimVariableNames <- function(openNc) {
  ncInfo <- file.inq.nc(openNc)
  varInfo <- lapply(seq_len(ncInfo$nvars)-1, var.inq.nc, ncfile = openNc)
  varNames <- sapply(varInfo, \(vi) vi$name)
  dimInfo <- lapply(seq_len(ncInfo$ndims)-1, dim.inq.nc, ncfile = openNc)
  dimNames <- sapply(dimInfo, \(vi) vi$name)
  return(setdiff(varNames, dimNames))
}

#' @export
ncGetVariableNames <- function(openNc) {
  ncInfo <- file.inq.nc(openNc)
  varInfo <- lapply(seq_len(ncInfo$nvars)-1, var.inq.nc, ncfile = openNc)
  varNames <- sapply(varInfo, \(vi) vi$name)
  return(varNames)
}

#' @export
ncGetDimensionNames <- function(openNc) {
  ncInfo <- file.inq.nc(openNc)
  dimInfo <- lapply(seq_len(ncInfo$ndims)-1, dim.inq.nc, ncfile = openNc)
  dimNames <- sapply(dimInfo, \(vi) vi$name)
  return(dimNames)
}

#' @export
ncGetDimensionIndex <- function(openNc, dimName) {
  ncInfo <- file.inq.nc(openNc)
  dimInfo <- lapply(seq_len(ncInfo$ndims)-1, dim.inq.nc, ncfile = openNc)
  dimNames <- sapply(dimInfo, \(vi) vi$name)
  idx <- which(dimName == dimNames)
  stopifnot(length(idx) == 1)
  return(idx)
}


openNc <- function(filePath) {
  if (!file.exists(filePath)) stop(filePath, " does not exist")
  open.nc(filePath)
}


varGetNc <- function(
  ncfile, variable, start=NA, count=NA,
  na.mode=NULL, collapse=TRUE, unpack=FALSE, rawchar=FALSE, fitnum=FALSE,
  cache_bytes=NA, cache_slots=NA, cache_preemption=NA) {
    if (is.null(na.mode)) {
      na.mode <- 4
    }
    var.get.nc(
      ncfile=ncfile, variable=variable, start=start, count=count,
      na.mode=na.mode, collapse=collapse, unpack=unpack, rawchar=rawchar, fitnum=fitnum,
      cache_bytes=cache_bytes, cache_slots=cache_slots, cache_preemption=cache_preemption
    )
}
