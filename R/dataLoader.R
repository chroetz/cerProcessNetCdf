loadData <- function(dataDescriptor) {
  subclass <- ConfigOpts::getClassAt(dataDescriptor, 2)
  switch(
    subclass,
    MultiFile = loadDataMultiFile(dataDescriptor),
    YearlyFiles = loadDataYearlyFiles(dataDescriptor),
    SingleFile = loadDataSingleFile(dataDescriptor),
    SingleFileTimeless = loadDataSingleFileTimeless(dataDescriptor),
    LabelFileTimeless = loadDataLabelFileTimeless(dataDescriptor),
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
  if (length(filePaths) == 0) {
    stop("No files found with DataDescriptor:", formatDataDesctiptor(dataDescriptor))
  }
  fileNames <- basename(filePaths)

  nc <- openNc(filePaths[1])
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
  startTime <- lapply(timeList, \(x) min(x)) |> unlist()
  endTime <- lapply(timeList, \(x) max(x)) |> unlist()
  names(timeList) <- filePaths
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
        startTime = startTime,
        endTime = endTime,
        year = yearList,
        fileName = fileNames,
        filePath = filePaths
      ) |>
        arrange(startTime) |>
        tidyr::unnest_longer(year)
    )

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
      descriptor = dataDescriptor,
      joinByLabel = FALSE,
      gridFormat,
      years = years,
      labels = labels,
      timeDimName,
      timeList,
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
  if (length(filePaths) == 0) {
    stop("No files found with DataDescriptor:", formatDataDesctiptor(dataDescriptor))
  }
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

  nc <- openNc(filePaths[1])
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
      joinByLabel = FALSE,
      gridFormat,
      years = unique(fileYears),
      labels = unique(fileLabels),
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


loadDataLabelFileTimeless <- function(dataDescriptor) {
  filePaths <- list.files(
    dataDescriptor$dirPath,
    pattern = dataDescriptor$pattern,
    recursive = dataDescriptor$recursive,
    full.names = TRUE)
  if (length(filePaths) == 0) {
    stop("No files found with DataDescriptor:", formatDataDesctiptor(dataDescriptor))
  }
  fileNames <- basename(filePaths)
  matchMatrix <- str_match(fileNames, dataDescriptor$pattern)
  if (ncol(matchMatrix) == 1) {
    fileLabels <- cerUtility::removeFileNameEnding(fileNames)
  } else {
    labelParts <- matchMatrix[,2:ncol(matchMatrix), drop=FALSE]
    fileLabels <- apply(labelParts, 1, paste, collapse="_")
  }

  nc <- openNc(filePaths[1])
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
      joinByLabel = TRUE,
      gridFormat,
      labels = unique(fileLabels),
      variableName,
      dimIds,
      dimNames,
      varDimIds,
      meta = tibble(
        label = fileLabels,
        fileName = fileNames,
        filePath = filePaths))
}


loadDataSingleFile <- function(dataDescriptor) {

  nc <- openNc(dataDescriptor$filePath)
  on.exit(close.nc(nc))

  dimNames <- ncGetDimensionNames(nc)
  timeDim <- ncLoadTimeDimension(nc)

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
    dim.inq.nc(nc, timeDim$name)$id)
  names(dimIds) <- c("lon", "lat", timeDim$name)
  dimNames <- c(
    dim.inq.nc(nc, 0)$name,
    dim.inq.nc(nc, 1)$name,
    dim.inq.nc(nc, 2)$name)

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
    joinByLabel = FALSE,
    descriptor = dataDescriptor,
    gridFormat,
    years = timeDim$years,
    labels,
    timeDimName = timeDim$name,
    dimIds,
    dimNames,
    varDimIds,
    meta = tidyr::expand_grid(label = labels, year = timeDim$years))
}


loadDataSingleFileTimeless <- function(dataDescriptor) {

  nc <- openNc(dataDescriptor$filePath)
  on.exit(close.nc(nc))

  dimNames <- ncGetDimensionNames(nc)

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
    dim.inq.nc(nc, "lat")$id)
  names(dimIds) <- c("lon", "lat")
  dimNames <- c(
    dim.inq.nc(nc, 0)$name,
    dim.inq.nc(nc, 1)$name)

  if (!"data" %in% names(.info)) .info$data <- list()
  .info$data[[dataDescriptor$name]] <- lst(
    joinByLabel = TRUE,
    descriptor = dataDescriptor,
    gridFormat,
    labels,
    dimIds,
    dimNames,
    varDimIds,
    meta = tidyr::expand_grid(label = labels))
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
    \(name) {
      lbls <- if (.info$data[[name]]$joinByLabel) labels else labels[[name]]
      env[[name]] <- getData(
        name,
        year,
        lbls,
        bbInfo)
    })
  return(invisible())
}


getData <- function(name, year, label = NULL, bbInfo = NULL) {

  dataInfo <- .info$data[[name]]

  if(!hasValue(label)) {
    label <- dataInfo$labels
  }
  if (length(label) > 1) { # Find best matching label
    labelIntersect <- intersect(label, dataInfo$labels)
    if (length(labelIntersect)==1) {
      label <- labelIntersect
    } else {
      labels <- stringr::str_split(label, "_") |> unlist()
      label <- intersect(labels, dataInfo$labels)
    }
  }
  if (length(label) != 1) {
    cat("\ndataInfo$labels:\n")
    dataInfo$labels |> print()
    cat("\nlabel:\n")
    label |> print()
    stop("label in getData is not a single value")
  }

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
    MultiFile = getDataMultiFile(dataInfo, year, label, bbInfoScaled),
    YearlyFiles = getDataYearlyFiles(dataInfo, year, label, bbInfoScaled),
    SingleFile = getDataSingleFile(dataInfo, year, label, bbInfoScaled),
    SingleFileTimeless = getDataSingleFileTimeless(dataInfo, year, label, bbInfoScaled),
    LabelFileTimeless = getDataLabelFileTimeless(dataInfo, year, label, bbInfoScaled),
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

  nc <- openNc(info$filePath)
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

  nc <- openNc(info$filePath)
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



getDataLabelFileTimeless <- function(dataInfo, year, label, bbInfo = NULL) {

  # NOTE: using filter() here seems a bit slow
  sel <- dataInfo$meta$label == label
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

  nc <- openNc(info$filePath)
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


getDataSingleTimeless <- function(dataInfo, year, label, bbInfo = NULL) {

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

  nc <- openNc(dataInfo$descriptor$filePath)

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

  nc <- openNc(dataInfo$descriptor$filePath)

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
  years <- Reduce(
    \(x, y) if (is.null(y)) x else intersect(x, y),
    yearsList)
  return(years)
}


getDataLabelsAndYearsAll <- function(yearsFilter) {
  nms <- names(.info$data)
  labelsAndYearsListAll <- lapply(nms, getDataLabelsAndYears, yearsFilter)
  labelsAndYearsList <- labelsAndYearsListAll
  names(labelsAndYearsList) <- nms
  isWithYear <- sapply(labelsAndYearsList, \(x) "year" %in% names(x))
  labelsAndYearsList <- labelsAndYearsList[isWithYear]
  nms <- nms[isWithYear]
  for (nm in nms) {
    names(labelsAndYearsList[[nm]])[names(labelsAndYearsList[[nm]]) == "label"] <- nm
  }
  labelsAndYears <- labelsAndYearsList |> first()
  for (nm in nms[-1]) {
    labelsAndYears <- inner_join(labelsAndYears, labelsAndYearsList[[nm]], join_by(year))
  }
  if (NROW(labelsAndYears) == 0) {
    cat("\ngetDataLabelsAndYearsAll(): labelsAndYearsListAll:\n")
    print(labelsAndYearsListAll)
    stop("No data found for the given years")
  }
  return(labelsAndYears)
}


getDataLabelsAndYears <- function(name, yearsFilter) {
  info <- .info$data[[name]]$meta
  hasTime <- "year" %in% names(info)
  if (hasTime) {
    labelsAndYears <-
      info |>
      select(.data$label, .data$year)
    if (hasValue(yearsFilter)) {
      labelsAndYears <- filter(labelsAndYears, .data$year %in% yearsFilter)
    }
  } else {
    labelsAndYears <-
      info |>
      select(.data$label)
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


formatDataDesctiptor <- function(dd) {
  sprintf(
    "path: %s, pattern: %s, recrusive: %s",
    dd$dirPath, dd$pattern, as.character(dd$recursive))
}
