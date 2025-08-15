getTheDataTimeAll <- function(lotIdx, lineCount, timeRange) {

  stopifnot(length(.info$data) == 1)
  dataInfo <- .info$data[[1]]
  label <- dataInfo$labels
  stopifnot(length(label) == 1)

  subclass <- ConfigOpts::getClassAt(dataInfo$descriptor, 2)
  data <- switch(
    subclass,
    MultiFile = getDataTimeAllMultiFile(dataInfo, label, lotIdx, lineCount, timeRange),
    #YearlyFiles = getDataTimeAllYearlyFiles(dataInfo, label, lotIdx, lineCount, timeRange), # TODO
    #SingleFile = getDataTimeAllSingleFile(dataInfo, label, lotIdx, lineCount, timeRange), # TODO
    stop("Unknown DataDescriptor subclass: ", subclass)
  )

  stopifnot(length(dim(data)) == 3)

  # TODO: switch lat order if necessary

  return(data)
}


getDataTimeAllMultiFile <- function(dataInfo, label, lotIdx, lineCount, timeRange = NULL) {

  sel <- dataInfo$meta$label == label
  if (!is.null(timeRange)) {
    yearFrom <- lubridate::year(timeRange[1])
    yearTo <- lubridate::year(timeRange[2])
    sel <- sel & dataInfo$meta$year >= yearFrom & dataInfo$meta$year <= yearTo
  }
  info <- dataInfo$meta[sel,]
  filePaths <- info$filePath |> unique()
  stopifnot(length(filePaths) > 0)

  dataList <- lapply(
    filePaths,
    getDataTimeAllMultiFileOneFile,
    dataInfo = dataInfo,
    label = label,
    lotIdx = lotIdx,
    lineCount = lineCount,
    timeRange = timeRange)

  data <- do.call(abind::abind, c(dataList, list(along = 3))) # time is the third dimension
  dimnames(data) <- dimnames(dataList[[1]]) # assumes no names for the values

  return(data)
}


getDataTimeAllMultiFileOneFile <- function(filePath, dataInfo, label, lotIdx, lineCount, timeRange = NULL) {

  times <- dataInfo$timeList[[filePath]]
  if (!is.null(timeRange)) {
    sel <- times >= timeRange[1] & times <= timeRange[2]
    timeIdx <- which(sel)
    stopifnot(all(abs(diff(timeIdx)) == 1))
  } else {
    timeIdx <- seq_along(times)
  }

  if (.info$idxDim == "lon") {
    lonLatTimeStart <- c(lotIdx, 1, min(timeIdx))
    lonLatTimeCount <- c(lineCount, NA, length(timeIdx))
  } else {
    lonLatTimeStart <- c(1, lotIdx, min(timeIdx))
    lonLatTimeCount <- c(NA, lineCount, length(timeIdx))
  }

  start <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeStart)
  count <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeCount)

  nc <- openNc(filePath)
  data <- varGetNc(
    nc,
    label,
    start = start,
    count = count,
    na.mode = dataInfo$descriptor$naMode,
    collapse = FALSE
  )
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
