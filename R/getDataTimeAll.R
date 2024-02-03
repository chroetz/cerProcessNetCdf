getDataTimeAll <- function(name, lotIdx, lineCount) {

  # TODO: maybe switching lon and lat makes it much faster??

  dataInfo <- .info$data[[name]]
  label <- dataInfo$labels
  stopifnot(length(label) == 1)

  subclass <- ConfigOpts::getClassAt(dataInfo$descriptor, 2)
  data <- switch(
    subclass,
    MultiFile = getDataTimeAllMultiFile(dataInfo, label, lotIdx, lineCount),
    YearlyFiles = getDataTimeAllYearlyFiles(dataInfo, label, lotIdx, lineCount), # TODO
    SingleFile = getDataTimeAllSingleFile(dataInfo, label, lotIdx, lineCount), # TODO
    stop("Unknown DataDescriptor subclass: ", subclass)
  )

  stopifnot(length(dim(data)) == 3)

  # TODO: switch lat order if necessary

  return(data)
}


getDataTimeAllMultiFile <- function(dataInfo, label, lotIdx, lineCount) {

  sel <- dataInfo$meta$label == label
  info <- dataInfo$meta[sel,]
  filePaths <- info$filePath |> unique()
  stopifnot(length(filePaths) > 0)

  dataList <- lapply(
    filePaths,
    getDataTimeAllMultiFileOneFile,
    dataInfo = dataInfo,
    label = label,
    lotIdx = lotIdx,
    lineCount = lineCount)

  data <- do.call(abind::abind, c(dataList, list(along = 3))) # time is the third dimension
  dimnames(data) <- dimnames(dataList[[1]]) # assumes no names for the values

  return(data)
}


getDataTimeAllMultiFileOneFile <- function(dataInfo, label, lotIdx, lineCount, filePath) {

  if (.info$idxDim == "lon") {
    lonLatTimeStart <- c(lotIdx, 1, 1)
    lonLatTimeCount <- c(lineCount, NA, NA)
  } else {
    lonLatTimeStart <- c(1, lotIdx, 1)
    lonLatTimeCount <- c(NA, lineCount, NA)
  }

  start <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeStart)
  count <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeCount)

  nc <- open.nc(filePath)
  data <- var.get.nc(
    nc,
    label,
    start = start,
    count = count,
    collapse = FALSE)
  close.nc(nc)

  if (dataInfo$descriptor$setNaToZero) {
    data <- ifelse(is.na(data), 0, data)
  }

  # Set correct dimnames of data
  dimNames <- dataInfo$dimNames[dataInfo$varDimIds+1]
  dimNameList <- list(NULL, NULL, NULL)
  names(dimNameList) <- dimNames
  dimnames(data) <- dimNameList

  return(data)
}
