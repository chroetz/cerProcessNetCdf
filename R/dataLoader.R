loadData <- function(name, dataDescriptor, targetVariable=NULL) {
  switch(
    dataDescriptor$type,
    YearlyFiles = loadDataYearlyFiles(name, dataDescriptor, targetVariable),
    SingleFile = loadDataSingleFile(name, dataDescriptor, targetVariable),
    stop("Unknown data descriptor type: ", dataDescriptor$type)
  )
  return(invisible())
}


loadDataYearlyFiles <- function(name, dataDescriptor, targetVariable = NULL) {
  fileNames <- list.files(dataDescriptor$dir, pattern=dataDescriptor$pattern)
  fileYears <- stringr::str_match(fileNames, dataDescriptor$pattern)[,2] |> as.integer()
  filePaths <- file.path(dataDescriptor$dir, fileNames)

  lonLatCheck <- checkLonLatInNcFile(filePaths[1])

  nc <- open.nc(filePaths[1])
  variableName <- ncGetNonDimVariableNames(nc)
  if (hasValue(targetVariable)) {
    variableName <- intersect(variableName, targetVariable)
  }
  stopifnot(length(variableName) == 1)
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
  .info$data[[name]] <- list(
    name = name,
    descriptor = dataDescriptor,
    lonIncreasing = lonLatCheck$lonIncreasing,
    latIncreasing = lonLatCheck$latIncreasing,
    years = fileYears,
    variableName = variableName,
    dimIds = dimIds,
    dimNames = dimNames,
    varDimIds = varDimIds,
    meta = tibble(
      year = fileYears,
      name = fileNames,
      path = filePaths))
}


loadDataSingleFile <- function(name, dataDescriptor, targetVariable = NULL) {

  lonLatCheck <- checkLonLatInNcFile(dataDescriptor$filePath)

  nc <- open.nc(dataDescriptor$filePath)

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
    stopifnot(stringr::str_detect(timeUnitDescription, pattern))
    startDayText <- stringr::str_match(timeUnitDescription, pattern)[,2]
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
  if (hasValue(targetVariable)) {
    variableName <- intersect(variableName, targetVariable)
  }
  stopifnot(length(variableName) == 1)

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

  close.nc(nc)

  stopifnot(max(abs(years - round(years))) < sqrt(.Machine$double.eps))

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[name]] <- list(
    name = name,
    descriptor = dataDescriptor,
    lonIncreasing = lonLatCheck$lonIncreasing,
    latIncreasing = lonLatCheck$latIncreasing,
    years = years,
    timeDimName = timeDimName,
    variableName = variableName,
    dimIds = dimIds,
    dimNames = dimNames,
    varDimIds = varDimIds,
    filePath = dataDescriptor$filePath)
}


checkLonLatInNcFile <- function(filePath) {
  nc <- open.nc(filePath)
  dimNames <- ncGetDimensionNames(nc)
  stopifnot(c("lon", "lat") %in% dimNames)
  lonValues <- var.get.nc(nc, "lon")
  latValues <- var.get.nc(nc, "lat")
  close.nc(nc)
  if (all(diff(lonValues) > 0)) lonIncreasing <- TRUE
  else if (all(diff(lonValues) < 0)) lonIncreasing <- FALSE
  else stop("Longitude values are neither increasing nor decreasing.")
  if (all(diff(latValues) > 0)) latIncreasing <- TRUE
  else if (all(diff(latValues) < 0)) latIncreasing <- FALSE
  else stop("Latitude values are neither increasing nor decreasing.")
  assertLonLat(
    if (lonIncreasing) lonValues else rev(lonValues),
    if (latIncreasing) latValues else rev(latValues))
  return(lst(lonIncreasing, latIncreasing))
}


getData <- function(name, year, setNaToZero = FALSE) {

  dataInfo <- .info$data[[name]]
  dataDescriptor <- dataInfo$descriptor
  data <- switch(
    dataDescriptor$type,
    YearlyFiles = getDataYearlyFiles(dataInfo, year, setNaToZero),
    SingleFile = getDataSingleFile(dataInfo, year, setNaToZero),
    stop("Unknown data descriptor type: ", dataDescriptor$type)
  )

  # Make sure that lon and lat are increasing
  if (!dataInfo$lonIncreasing == .info$grid$lonIncreasing) {
    data <- reverseArrayDim(data, which(dimnames(data) |> names() == "lon"))
  }
  if (!dataInfo$latIncreasing == .info$grid$latIncreasing) {
    data <- reverseArrayDim(data, which(dimnames(data) |> names() == "lat"))
  }
  if (!all(data |> dimnames() |> names() == .info$grid$orderedNames)) {
    data <- t.default(data)
  }

  return(data)
}


getDataYearlyFiles <- function(dataInfo, year, setNaToZero = FALSE) {

  fileInfo <- dataInfo$meta |> filter(.data$year == .env$year)
  stopifnot(nrow(fileInfo) == 1)

  nc <- open.nc(fileInfo$path)
  data <- var.get.nc(nc, dataInfo$variableName)
  close.nc(nc)

  if (setNaToZero) {
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


getDataSingleFile <- function(dataInfo, year, setNaToZero = FALSE) {

  timeIdx <- which(dataInfo$years == year)

  lonLatTimeStart <- c(1, 1, timeIdx)
  lonLatTimeCount <- c(NA, NA, 1)
  start <- permuteDimIds(dataInfo, lonLatTimeStart)
  count <- permuteDimIds(dataInfo, lonLatTimeCount)

  nc <- open.nc(dataInfo$filePath)

  data <- var.get.nc(
    nc,
    dataInfo$variableName,
    start = start,
    count = count
  )
  close.nc(nc)

  if (setNaToZero) {
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

