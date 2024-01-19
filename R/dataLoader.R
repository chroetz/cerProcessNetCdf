loadData <- function(dataDescriptor) {
  subclass <- ConfigOpts::getClassAt(dataDescriptor, 2)
  switch(
    subclass,
    YearlyFiles = loadDataYearlyFiles(dataDescriptor),
    SingleFile = loadDataSingleFile(dataDescriptor),
    stop("Unknown DataDescriptor subclass: ", subclass)
  )

  if (hasValue(.info$boundingBoxes)) {
    # TODO: currently bounding box is only supported when grid formats are identical
    stopifnot(identical(.info$targetFormat, .info$data[[dataDescriptor$name]]$gridFormat))
  }

  return(invisible())
}


loadDataYearlyFiles <- function(dataDescriptor) {
  fileNames <- list.files(dataDescriptor$dir, pattern=dataDescriptor$pattern)
  fileYears <- str_match(fileNames, dataDescriptor$pattern)[,2] |> as.integer()
  filePaths <- file.path(dataDescriptor$dir, fileNames)

  nc <- open.nc(filePaths[1])
  variableName <- ncGetNonDimVariableNames(nc)
  if (hasValueString(dataDescriptor$dataVariableName)) {
    variableName <- intersect(variableName, dataDescriptor$dataVariableName)
  }
  stopifnot(length(variableName) == 1)
  gridFormat <- getNativeGridFormatFromNc(nc, variableName)
  cat(
    "Grid format of variable", variableName, ":",
    format(gridFormat),
    "\n")
  varInfo <- var.inq.nc(nc, variableName)
  varDimIds <- varInfo$dimids
  dimIds <- c(
    dim.inq.nc(nc, "lon")$id,
    dim.inq.nc(nc, "lat")$id)
  names(dimIds) <- c("lon", "lat")
  dimNames <- c(
    dim.inq.nc(nc, 0)$name,
    dim.inq.nc(nc, 1)$name)
  close.nc(nc)

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
      descriptor = dataDescriptor,
      gridFormat,
      years = fileYears,
      variableName,
      dimIds,
      dimNames,
      varDimIds,
      meta = tibble(
        year = fileYears,
        name = fileNames,
        path = filePaths))
}


loadDataSingleFile <- function(dataDescriptor) {

  nc <- open.nc(dataDescriptor$filePath)
  on.exit(close.nc(nc))

  dimNames <- ncGetDimensionNames(nc)
  timeDimName <- setdiff(dimNames, c("lon", "lat"))
  stopifnot(length(timeDimName) == 1)
  timeValues <- var.get.nc(nc, timeDimName)
  timeVarInfo <- var.inq.nc(nc, timeDimName)

  attNames <- sapply(
    seq_len(timeVarInfo$natts)-1,
    \(i) att.inq.nc(nc, timeDimName, i)$name)
  if ("units" %in% attNames) {
    timeUnitDescription <- att.get.nc(nc, timeDimName, "units")
    pattern <- "^days since ([\\d-]+)( \\d{2}:\\d{2}:(\\d{2})?)?"
    stopifnot(str_detect(timeUnitDescription, pattern))
    startDayText <- str_match(timeUnitDescription, pattern)[,2]
    startDate <- as.Date(startDayText)
    startYear <- lubridate::year(startDate)
    years <- timeValues/365 + startYear
    formattedStartDate <- format(startDate, "%B %d, %Y")
    cat(
      "Assume that time values are days since year", startYear, "(", formattedStartDate, ").\n")
  } else {
    years <- timeValues
    cat("Assume that time values are years.\n")
  }

  variableName <- ncGetNonDimVariableNames(nc)
  if (hasValueString(dataDescriptor$dataVariableName)) {
    variableName <- intersect(variableName, dataDescriptor$dataVariableName)
  }
  stopifnot(length(variableName) == 1)
  gridFormat <- getNativeGridFormatFromNc(nc, variableName)
  cat(
    "Grid format of variable", variableName, ":",
    format(gridFormat),
    "\n")

  varInfo <- var.inq.nc(nc, variableName)
  varDimIds <- varInfo$dimids
  dimIds <- c(
    dim.inq.nc(nc, "lon")$id,
    dim.inq.nc(nc, "lat")$id,
    dim.inq.nc(nc, timeDimName)$id)
  names(dimIds) <- c("lon", "lat", timeDimName)
  dimNames <- c(
    dim.inq.nc(nc, 0)$name,
    dim.inq.nc(nc, 1)$name,
    dim.inq.nc(nc, 2)$name)

  stopifnot(max(abs(years - round(years))) < sqrt(.Machine$double.eps))

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
    descriptor = dataDescriptor,
    years,
    timeDimName,
    variableName,
    dimIds,
    dimNames,
    varDimIds,
    gridFormat)
}


getDataAll <- function(year, bbInfo = NULL) {
  data <- lapply(
    names(.info$data),
    \(name) getData(name, year, bbInfo))
  names(data) <- names(.info$data)
  return(data)
}


getData <- function(name, year, bbInfo = NULL) {

  dataInfo <- .info$data[[name]]
  subclass <- ConfigOpts::getClassAt(dataInfo$descriptor, 2)
  data <- switch(
    subclass,
    YearlyFiles = getDataYearlyFiles(dataInfo, year, bbInfo),
    SingleFile = getDataSingleFile(dataInfo, year, bbInfo),
    stop("Unknown DataDescriptor subclass: ", subclass)
  )

  # Make sure that lon and lat are correctly ordered
  if (dataInfo$gridFormat$lonIncreasing != .info$grid$lonIncreasing) {
    data <- reverseArrayDim(data, which(dimnames(data) |> names() == "lon"))
  }
  if (dataInfo$gridFormat$latIncreasing != .info$grid$latIncreasing) {
    data <- reverseArrayDim(data, which(dimnames(data) |> names() == "lat"))
  }
  if (!all(data |> dimnames() |> names() == .info$grid$orderedNames)) {
    data <- t.default(data)
  }

  return(data)
}


getDataYearlyFiles <- function(dataInfo, year, bbInfo = NULL) {

  fileInfo <- dataInfo$meta |> filter(.data$year == .env$year)
  stopifnot(nrow(fileInfo) == 1)

  if (hasValue(bbInfo)) {
    # TODO: convert bbox format to data format
    lonLatStart <- c(
      bbInfo$min_lon,
      bbInfo$min_lat)
    lonLatCount <- c(
      bbInfo$max_lon - bbInfo$min_lon + 1,
      bbInfo$max_lat - bbInfo$min_lat + 1)
  } else {
    lonLatStart <- c(1, 1)
    lonLatCount <- c(NA, NA)
  }
  start <- permuteDimIdsLonLat(dataInfo, lonLatStart)
  count <- permuteDimIdsLonLat(dataInfo, lonLatCount)

  nc <- open.nc(fileInfo$path)
  data <- var.get.nc(
    nc,
    dataInfo$variableName,
    start = start,
    count = count)
  close.nc(nc)

  if (dataInfo$descriptor$setNaToZero) {
    data <- ifelse(is.na(data), 0, data)
  }

  # Set correct dimnames of data
  dimNames <- dataInfo$dimNames[dataInfo$varDimIds+1]
  dimNames <- dimNames[dimNames %in% c("lon", "lat")]
  dimNameList <- list(NULL, NULL)
  names(dimNameList) <- dimNames
  dimnames(data) <- dimNameList

  return(data)
}


getDataSingleFile <- function(dataInfo, year, bbInfo = NULL) {

  timeIdx <- which(dataInfo$years == year)

  if (hasValue(bbInfo)) {
    # TODO: convert bbox format to data format
    lonLatTimeStart <- c(
      bbInfo$min_lon,
      bbInfo$min_lat,
      timeIdx)
    lonLatTimeCount <- c(
      bbInfo$max_lon - bbInfo$min_lon + 1,
      bbInfo$max_lat - bbInfo$min_lat + 1,
      1)
  } else {
    lonLatTimeStart <- c(1, 1, timeIdx)
    lonLatTimeCount <- c(NA, NA, 1)
  }
  start <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeStart)
  count <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeCount)

  nc <- open.nc(dataInfo$descriptor$filePath)

  data <- var.get.nc(
    nc,
    dataInfo$variableName,
    start = start,
    count = count
  )
  close.nc(nc)

  if (dataInfo$descriptor$setNaToZero) {
    data <- ifelse(is.na(data), 0, data)
  }

  # Set correct dimnames of data
  dimNames <- dataInfo$dimNames[dataInfo$varDimIds+1]
  dimNames <- dimNames[dimNames %in% c("lon", "lat")]
  dimNameList <- list(NULL, NULL)
  names(dimNameList) <- dimNames
  dimnames(data) <- dimNameList

  return(data)
}


permuteDimIdsLonLatTime <- function(dataInfo, lonLatTimeVec) {
  if (is.null(names(lonLatTimeVec))) {
    names(lonLatTimeVec) <- c("lon", "lat", dataInfo$timeDimName)
  }
  permuter <- dataInfo$dimNames[dataInfo$varDimIds+1]
  permutedVec <- lonLatTimeVec[permuter]
  return(permutedVec)
}


permuteDimIdsLonLat <- function(dataInfo, lonLatVec) {
  if (is.null(names(lonLatVec))) {
    names(lonLatVec) <- c("lon", "lat")
  }
  permuter <- dataInfo$dimNames[dataInfo$varDimIds+1]
  permutedVec <- lonLatVec[permuter]
  return(permutedVec)
}


permuteDimIds <- function(dataInfo, lonLatTimeVec) {
  if (is.null(names(lonLatTimeVec))) {
    names(lonLatTimeVec) <- c("lon", "lat", dataInfo$timeDimName)
  }
  permuter <- dataInfo$dimNames[dataInfo$varDimIds+1]
  permutedVec <- lonLatTimeVec[permuter]
  return(permutedVec)
}


getDataYears <- function(name) {
  dataInfo <- .info$data[[name]]
  return(dataInfo$years)
}


getDataYearsAll <- function() {
  yearsList <- lapply(names(.info$data), getDataYears)
  years <- Reduce(intersect, yearsList)
  return(years)
}

