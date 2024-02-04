rollTimeApply <- function(
  targetFormat = NULL,
  variableDataDescriptor,
  applyFunctionList,
  outFilePath,
  padding = c(3,3),
  fill = NA_real_,
  idxDim = c("lon", "lat"),
  lineCount = 1,
  timeRange = NULL,
  nBatches = 1,
  batchIndex = 1
) {

  idxDim <- match.arg(idxDim)

  clearInfo()

  loadData(variableDataDescriptor)

  if (is.null(targetFormat)) targetFormat <- .info$data[[1]]$gridFormat

  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(argNames, \(nm) assign(nm, env[[nm]], .info))

  initializeGrid(targetFormat)

  if (file.exists(.info$outFilePath)) {
    cat(.info$outFilePath, "already exists. Deleting.\n")
    file.remove(.info$outFilePath)
  }

  lotMax <- if (idxDim == "lon") .info$grid$nLon else .info$grid$nLat
  lonIdices <- seq(1, lotMax, by = lineCount)
  lineCounts <- rep(lineCount, length(lonIdices))
  lineCounts[length(lonIdices)] <- lotMax - lonIdices[length(lonIdices)] + 1
  allBatchElements <- lapply(seq_along(lonIdices), \(i) list(lotIdx = lonIdices[i], lineCount = lineCounts[i]))
  batch <- cerUtility::splitAndGetOneBatch(
    paste0(idxDim, "-index"),
    allBatchElements,
    nBatches,
    batchIndex)

  cat("Start main loop.\n")
  for (batchElement in batch) {
    processRollTimeApply(batchElement$lotIdx, batchElement$lineCount, timeRange)
  }
  cat("End main loop.\n")
}


processRollTimeApply <- function(lotIdx, lineCount, timeRange) {

  pt <- proc.time()
  cat("Processing index", lotIdx, "...\n")

  pt1 <- proc.time()
  data <- getTheDataTimeAll(lotIdx, lineCount, timeRange)
  cat("\tgetDataAllTimeAll duration:", (proc.time()-pt1)[3], "s\n")

  pt2 <- proc.time()
  valueList <- lapply(
    .info$applyFunctionList,
    \(fun) {
      res <- array(NA_real_, dim = dim(data), dimnames = dimnames(data))
      for (iLon in seq_len(dim(data)[1])) # assumes order: lon lat time
        for (iLat in seq_len(dim(data)[2])) {
          d <- data[iLon, iLat, ]
          res[iLon, iLat, ] <- rollApply(fun, d, .info$padding, .info$fill)
        }
      return(res)
    })
  names(valueList) <- names(.info$applyFunctionList)
  cat("\tapplyFunctionList duration:", (proc.time()-pt2)[3], "s\n")

  outFilePath <- paste0(
    cerUtility::removeFileNameEnding(.info$outFilePath),
    "_", lotIdx, ".nc")

  if (.info$idxDim == "lon") {
    dimList <- list(
      lon = .info$grid$lonValues[lotIdx:(lotIdx + lineCount - 1)],
      lat = .info$grid$latValues)
  } else {
    dimList <- list(
      lon = .info$grid$lonValues,
      lat = .info$grid$latValues[lotIdx:(lotIdx + lineCount - 1)])
  }
  allTimes <- .info$data[[1]]$timeList |> unlist()
  times <- allTimes[allTimes >= .info$timeRange[1] & allTimes <= .info$timeRange[2]]
  dimList[[.info$data[[1]]$timeDimName]] <- as.integer(times)

  pt3 <- proc.time()
  saveNetCdf(
    outFilePath,
    dimList = dimList,
    valueList = valueList)
  cat("\tSaving to", outFilePath, "duration:", (proc.time()-pt3)[3], "s\n")
  cat("All done after", (proc.time()-pt)[3], "s\n")
}


rollApply <- function(fun, x, padding, fill = NA_real_) {
  n <- length(x)
  if (is.character(fun)) {
    width <- sum(padding) + 1
    resApply <- switch(
      fun,
      mean = zoo::rollmean(x, width),
      max = zoo::rollmax(x, width),
      median = zoo::rollmedian(x, width),
      sum = zoo::rollsum(x, width),
      stop("Unknown function string ", fun, "."))
  } else if (is.function(fun)) {
    resApply <- vapply(
      (1+padding[1]):(n-padding[2]),
      \(i) fun(x[(i-padding[1]):(i+padding[2])]),
      numeric(1))
  } else {
    stop("fun must be function or character.")
  }

  if (is.numeric(fill)) {
    res <- c(rep(fill, padding[1]), resApply, rep(fill, padding[2]))
  } else if (is.character(fill)) {
    res <- switch(
      fill,
      const = c(rep(first(resApply),padding[1]), resApply, rep(last(resApply), padding[2])),
      stop("Unknown fill string ", fill, "."))
  } else {
    stop("fill must be numeric or character.")
  }
  return(res)
}
