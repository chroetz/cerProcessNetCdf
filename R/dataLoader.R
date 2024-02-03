loadData <- function(dataDescriptor) {
  subclass <- ConfigOpts::getClassAt(dataDescriptor, 2)
  switch(
    subclass,
    MultiFile = loadDataMultiFile(dataDescriptor),
    YearlyFiles = loadDataYearlyFiles(dataDescriptor),
    SingleFile = loadDataSingleFile(dataDescriptor),
    stop("Unknown DataDescriptor subclass: ", subclass)
  )

  return(invisible())
}


loadDataMultiFile <- function(dataDescriptor) {

  filePaths <- list.files(
    dataDescriptor$dirPath,
    pattern = dataDescriptor$pattern,
    recursive = dataDescriptor$recursive,
    full.names = TRUE)
  fileNames <- basename(filePaths)

  nc <- open.nc(filePaths[1])
  labels <- ncGetNonDimVariableNames(nc)
  if (hasValueString(dataDescriptor$dataVariableNames)) {
    labels <- intersect(labels, dataDescriptor$dataVariableNames)
  }
  gridFormat <- getNativeGridFormatFromNc(nc, labels[1])
  cat(
    "Grid format of variable", labels[1], ":",
    format(gridFormat),
    "\n")
  if (length(labels) > 1) {
    cat("WARNING: Found more than one variable in file. Assume they all have the same grid format\n")
  }

  dimNames <- ncGetDimensionNames(nc)
  timeDimName <- setdiff(dimNames, c("lon", "lat"))
  stopifnot(length(timeDimName) == 1)
  varInfo <- var.inq.nc(nc, labels[1])
  varDimIds <- varInfo$dimids
  dimIds <- c(
    dim.inq.nc(nc, "lon")$id,
    dim.inq.nc(nc, "lat")$id,
    dim.inq.nc(nc, timeDimName)$id)
  names(dimIds) <- c("lon", "lat", timeDimName)
  close.nc(nc)

  timeList <- lapply(filePaths, getFileTimes, timeDimName = timeDimName)
  yearList <- lapply(timeList, \(x) lubridate::year(x) |> unique())
  times <- unlist(timeList)
  years <- unlist(yearList)
  stopifnot(length(times) == length(unique(times)))
  stopifnot(length(years) == length(unique(years)))

  meta <-
    cross_join(
      tibble(
        label = labels),
      tibble(
        year = yearList,
        fileName = fileNames,
        filePath = filePaths
      ) |>
      tidyr::unnest_longer(year)
    )

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
      descriptor = dataDescriptor,
      gridFormat,
      years = years,
      label = unique(fileLabels),
      variableName,
      timeDimName,
      dimIds,
      dimNames,
      varDimIds,
      meta)
}



loadDataYearlyFiles <- function(dataDescriptor) {
  filePaths <- list.files(
    dataDescriptor$dirPath,
    pattern = dataDescriptor$pattern,
    recursive = dataDescriptor$recursive,
    full.names = TRUE)
  fileNames <- basename(filePaths)
  matchMatrix <- str_match(fileNames, dataDescriptor$pattern)
  # Assume that the last capture group is the year
  fileYears <- matchMatrix[,ncol(matchMatrix)] |> as.integer()
  if (ncol(matchMatrix) == 2) {
    fileLabels <- cerUtility::longestCommonPrefix(fileNames)
  } else {
    labelParts <- matchMatrix[,2:(ncol(matchMatrix)-1), drop=FALSE]
    fileLabels <- apply(labelParts, 1, paste, collapse="_")
  }

  nc <- open.nc(filePaths[1])
  variableName <- ncGetNonDimVariableNames(nc)
  if (hasValueString(dataDescriptor$dataVariableName)) {
    variableName <- intersect(variableName, dataDescriptor$dataVariableName)
  }
  stopifnot(length(variableName) == 1)
  gridFormat <- getNativeGridFormatFromNc(nc, variableName)
  ncInfo <- file.inq.nc(nc)
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
  dimNames <- sapply(seq_len(ncInfo$ndims)-1, \(dimid) dim.inq.nc(nc, dimid)$name)
  close.nc(nc)

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
      descriptor = dataDescriptor,
      gridFormat,
      years = unique(fileYears),
      label = unique(fileLabels),
      variableName,
      dimIds,
      dimNames,
      varDimIds,
      meta = tibble(
        label = fileLabels,
        year = fileYears,
        fileName = fileNames,
        filePath = filePaths))
}


loadDataSingleFile <- function(dataDescriptor) {

  nc <- open.nc(dataDescriptor$filePath)
  on.exit(close.nc(nc))

  dimNames <- ncGetDimensionNames(nc)
  timeDimName <- setdiff(dimNames, c("lon", "lat"))
  stopifnot(length(timeDimName) == 1)
  timeValues <- var.get.nc(nc, timeDimName) |> as.vector()
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

  labels <- ncGetNonDimVariableNames(nc)
  if (hasValueString(dataDescriptor$dataVariableNames)) {
    labels <- intersect(labels, dataDescriptor$dataVariableNames)
  }
  gridFormat <- getNativeGridFormatFromNc(nc, labels[1])
  cat(
    "Grid format of variable", labels[1], ":",
    format(gridFormat),
    "\n")
  if (length(labels) > 1) {
    cat("WARNING: Found more than one variable in file. Assume they all have the same grid format\n")
  }

  varInfo <- var.inq.nc(nc, labels[1])
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
    gridFormat,
    years,
    labels,
    timeDimName,
    dimIds,
    dimNames,
    varDimIds,
    meta = tidyr::expand_grid(label = labels, year = years))
}


getDataAll <- function(year, label = NULL, bbInfo = NULL) {
  data <- lapply(
    names(.info$data),
    \(name) getData(name, year, label, bbInfo))
  names(data) <- names(.info$data)
  return(data)
}

assignDataAll <- function(year, labels, env, bbInfo = NULL) {
  lapply(
    names(.info$data),
    \(name) env[[name]] <- getData(name, year, labels[[name]], bbInfo))
  return(invisible())
}


getData <- function(name, year, label = NULL, bbInfo = NULL) {

  dataInfo <- .info$data[[name]]

  if(!hasValue(label)) {
    label <- dataInfo$labels
  }
  stopifnot(length(label) == 1)

  if (hasValue(bbInfo)) {

    bbInfo <- convertBoundingBoxLonLatIncrDecr(bbInfo, .info$boundingBoxFormat,  dataInfo$gridFormat)

    bbInfoScaled <- list(
      min_lon = pmax(1, floor(bbInfo$min_lon / dataInfo$descriptor$blowUpLon)),
      max_lon = ceiling(bbInfo$max_lon / dataInfo$descriptor$blowUpLon),
      min_lat = pmax(1, floor(bbInfo$min_lat / dataInfo$descriptor$blowUpLat)),
      max_lat = ceiling(bbInfo$max_lat / dataInfo$descriptor$blowUpLat))
  } else {
    bbInfoScaled <- NULL
  }

  subclass <- ConfigOpts::getClassAt(dataInfo$descriptor, 2)
  data <- switch(
    subclass,
    YearlyFiles = getDataYearlyFiles(dataInfo, year, label, bbInfoScaled),
    SingleFile = getDataSingleFile(dataInfo, year, label, bbInfoScaled),
    stop("Unknown DataDescriptor subclass: ", subclass)
  )

  data <- ensureGridFormat(data, .info$grid, dataInfo$gridFormat)

  data <- blowUp(
    data,
    dataInfo$descriptor$blowUpLon,
    dataInfo$descriptor$blowUpLat,
    bbInfo, bbInfoScaled)

  return(data)
}


getDataMultiFile <- function(dataInfo, year, label, bbInfo = NULL) {

  # TODO

  # NOTE: using filter() here seems a bit slow
  sel <- dataInfo$meta$year == year & dataInfo$meta$label == label
  info <- dataInfo$meta[sel,]
  stopifnot(nrow(info) == 1)

  timeIndices <- getTimeIndicesYear(
    info$filePath,
    dataInfo$timeDimName,
    year)
  stopifnot(all(abs(diff(timeIndices)) == 1))


  if (hasValue(bbInfo)) {
    # TODO: convert bbox format to data format
    lonLatTimeStart <- c(
      bbInfo$min_lon,
      bbInfo$min_lat,
      min(timeIndices))
    lonLatTimeCount <- c(
      bbInfo$max_lon - bbInfo$min_lon + 1,
      bbInfo$max_lat - bbInfo$min_lat + 1,
      length(timeIndices))
  } else {
    lonLatTimeStart <- c(1, 1, min(timeIndices))
    lonLatTimeCount <- c(NA, NA, length(timeIndices))
  }
  start <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeStart)
  count <- permuteDimIdsLonLatTime(dataInfo, lonLatTimeCount)

  nc <- open.nc(info$filePath)
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


getDataYearlyFiles <- function(dataInfo, year, label, bbInfo = NULL) {

  # NOTE: using filter() here seems a bit slow
  sel <- dataInfo$meta$year == year & dataInfo$meta$label == label
  info <- dataInfo$meta[sel,]
  stopifnot(nrow(info) == 1)

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

  nc <- open.nc(info$filePath)
  data <- var.get.nc(
    nc,
    dataInfo$variableName,
    start = start,
    count = count,
    collapse = FALSE)
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


getDataSingleFile <- function(dataInfo, year, label, bbInfo = NULL) {

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
    label,
    start = start,
    count = count,
    collapse = FALSE
  )
  close.nc(nc)

  if (dataInfo$descriptor$setNaToZero) {
    data <- ifelse(is.na(data), 0, data)
  }

  # Set correct dimnames of data
  dimNames <- dataInfo$dimNames[dataInfo$varDimIds+1]
  dimSel <- dimNames %in% c("lon", "lat")
  dimNames <- dimNames[dimSel]
  dimNameList <- list(NULL, NULL)
  names(dimNameList) <- dimNames
  dim(data) <- dim(data)[dimSel]
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


getDataLabelsAndYearsAll <- function(yearsFilter) {
  nms <- names(.info$data)
  labelsAndYearsList <- lapply(nms, getDataLabelsAndYears, yearsFilter)
  names(labelsAndYearsList) <- nms
  for (nm in nms) {
    names(labelsAndYearsList[[nm]])[names(labelsAndYearsList[[nm]]) == "label"] <- nm
  }
  labelsAndYears <- labelsAndYearsList |> first()
  for (nm in nms[-1]) {
    labelsAndYears <- inner_join(labelsAndYears, labelsAndYearsList[[nm]], join_by(year))
  }
  return(labelsAndYears)
}


getDataLabelsAndYears <- function(name, yearsFilter) {
  labelsAndYears <-
    .info$data[[name]]$meta |>
    select(.data$label, .data$year)
  if (hasValue(yearsFilter)) {
    labelsAndYears <- filter(labelsAndYears, .data$year %in% yearsFilter)
  }
  return(labelsAndYears)
}


blowUp <- function(x, blowUpLon, blowUpLat, bbInfo, bbInfoScaled) {
  if (blowUpLon == 1 && blowUpLat == 1) return(x)
  stopifnot(is.matrix(x))
  dimNamesNames <- names(dimnames(x))
  stopifnot(length(dimNamesNames) == 2)

  lonOffset <- bbInfo$min_lon - ((bbInfoScaled$min_lon-1) * blowUpLon + 1)
  lonLength <- bbInfo$max_lon - bbInfo$min_lon + 1
  lonIdx <- lonOffset + seq_len(lonLength)
  latOffset <- bbInfo$min_lat - ((bbInfoScaled$min_lat-1) * blowUpLat + 1)
  latLength <- bbInfo$max_lat - bbInfo$min_lat + 1
  latIdx <- latOffset + seq_len(latLength)

  if (all(dimNamesNames == c("lon", "lat"))) {
    xRowBlown <- x[rep(seq_len(nrow(x)), each = blowUpLon), , drop=FALSE]
    xBlown <- xRowBlown[, rep(seq_len(ncol(x)), each = blowUpLat), drop=FALSE]
    y <- xBlown[lonIdx, latIdx, drop=FALSE]
  } else if (all(dimNamesNames == c("lat", "lon"))) {
    xRowBlown <- x[rep(seq_len(nrow(x)), each = blowUpLat), , drop=FALSE]
    xBlown <- xRowBlown[, rep(seq_len(ncol(x)), each = blowUpLon), drop=FALSE]
    y <- xBlown[latIdx, lonIdx, drop=FALSE]
  } else {
    stop("Unknown dimnames: ", dimNamesNames)
  }

  return(y)
}


