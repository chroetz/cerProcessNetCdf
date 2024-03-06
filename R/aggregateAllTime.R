#' @export
aggregateAllTime <- function(
  targetFormat = NULL,
  variableDataDescriptor,
  applyFunctionList,
  outFilePath,
  idxDim = c("lon", "lat"),
  lineCount = 1,
  timeRange = NULL,
  deflate = 9,
  nBatches = 1,
  batchIndex = 1
) {

  if (!is.null(timeRange)) timeRange <- lubridate::as_date(timeRange)

  idxDim <- match.arg(idxDim)

  applyFunctionListNames <- names(applyFunctionList)
  applyFunctionList <- lapply(applyFunctionList, resolveApplyFunction)
  names(applyFunctionList) <- applyFunctionListNames

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
  allBatchElements <- lapply(
    seq_along(lonIdices),
    \(i) list(lotIdx = lonIdices[i], lineCount = lineCounts[i]))
  batch <- cerUtility::splitAndGetOneBatch(
    paste0(idxDim, "-index"),
    allBatchElements,
    nBatches,
    batchIndex)

  cat("Start main loop.\n")
  for (batchElement in batch) {
    processAggregateAllTime(batchElement$lotIdx, batchElement$lineCount, timeRange)
  }
  cat("End main loop.\n")
}


processAggregateAllTime <- function(lotIdx, lineCount, timeRange) {

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
      # TODO: assumes order: lon lat time
      res <- array(NA_real_, dim = dim(data)[1:2], dimnames = dimnames(data)[1:2])
      for (iLon in seq_len(dim(data)[1]))
        for (iLat in seq_len(dim(data)[2])) {
          d <- data[iLon, iLat, ]
          res[iLon, iLat] <- fun(d)
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
