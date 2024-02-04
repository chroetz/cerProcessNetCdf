#' @export
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
  deflate = 9,
  nBatches = 1,
  batchIndex = 1
) {

  idxDim <- match.arg(idxDim)

  applyFunctionListNames <- names(applyFunctionList)
  applyFunctionList <- lapply(applyFunctionList, resolveAppltFunction)
  names(applyFunctionList) <- applyFunctionListNames

  if (!length(fill) == 0) fill <- NA_real_
  if (is.character(fill) && !fill %in% c("const")) fill <- eval(parse(text = fill))

  if (!is.null(timeRange)) timeRange <- lubridate::as_datetime(timeRange)

  clearInfo()
  loadData(variableDataDescriptor)
  if (is.null(targetFormat)) targetFormat <- .info$data[[1]]$gridFormat

  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(argNames, \(nm) assign(nm, env[[nm]], .info))

  initializeGrid(targetFormat)

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
  cat("\tGetting data ...")
  data <- getTheDataTimeAll(lotIdx, lineCount, timeRange)
  cat("duration:", (proc.time()-pt1)[3], "s\n")

  pt2 <- proc.time()
  cat("\tApplying functions ...")
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
  cat("duration:", (proc.time()-pt2)[3], "s\n")

  outFilePath <- paste0(
    cerUtility::removeFileNameEnding(.info$outFilePath),
    "_", lotIdx, ".nc")
  dirPath <- dirname(outFilePath)
  if (!dir.exists(dirPath)) dir.create(dirPath, recursive = TRUE)

  if (.info$idxDim == "lon") {
    dimList <- list(
      lon = .info$grid$lonValues[lotIdx:(lotIdx + lineCount - 1)],
      lat = .info$grid$latValues)
  } else {
    dimList <- list(
      lon = .info$grid$lonValues,
      lat = .info$grid$latValues[lotIdx:(lotIdx + lineCount - 1)])
  }

  timeValues <- .info$data[[1]]$timeValuesList |> unlist()
  if (!is.null(.info$timeRange)) {
    times <- .info$data[[1]]$timeList |> unlist()
    timeValues <- timeValues[times >= .info$timeRange[1] & times <= .info$timeRange[2]]
  }
  stopifnot(all(order(timeValues) == seq_along(timeValues))) # is sorted
  dimList[[.info$data[[1]]$timeDimName]] <- timeValues

  pt3 <- proc.time()
  cat("\tSaving to", outFilePath, "...")
  saveNetCdf(
    outFilePath,
    dimList = dimList,
    valueList = valueList,
    deflate = .info$deflate)
  cat("duration:", (proc.time()-pt3)[3], "s\n")
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


resolveAppltFunction <- function(fun) {
  if (is.function(fun)) return(fun)
  if (is.character(fun)) {
    if (fun %in% c("mean", "max", "median", "sum")) return(fun)
    return(eval(parse(text = fun)))
  } else {
    stop("fun must be function or character.")
  }
}
