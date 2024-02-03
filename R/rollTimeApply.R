rollTimeApply <- function(
  targetFormat = NULL,
  variableDataDescriptorList,
  applyExpressionList,
  outFilePath,
  idxDim = c("lon", "lat"),
  lineCount = 1,
  timeRange = NULL, #TODO
  nBatches = 1,
  batchIndex = 1
) {

  idxDim <- match.arg(idxDim)

  clearInfo()

  for (variableDataDescriptor in variableDataDescriptorList$list) {
    loadData(variableDataDescriptor)
  }

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
    processRollTimeApply(batchElement$lotIdx, batchElement$lineCount)
  }
  cat("End main loop.\n")
}


processRollTimeApply <- function(lotIdx, lineCount) {

  pt <- proc.time()
  cat("Processing index", lotIdx, "...\n")

  pt1 <- proc.time()
  dataList <- getDataAllTimeAll(lotIdx, lineCount)
  cat("\tgetDataAllTimeAll duration:", (proc.time()-pt1)[3], "s\n")

  dim3 <- dim(dataList[[1]])
  dimNames <- dimnames(dataList[[1]])

  pt2 <- proc.time()
  valueList <- lapply(
    .info$applyExpressionList,
    \(expr) {
      res <- array(NA_real_, dim = dim3, dimnames = dimNames)
      for (iLon in seq_len(dim3[1])) for (iLat in seq_len(dim3[2])) { # assume lon lat time order
        d <- lapply(dataList, \(x) x[iLon, iLat, ])
        names(d) <- names(dataList)
        res[iLon, iLat, ] <- rollEval(expr, d)
      }
      return(res)
    })
  names(valueList) <- names(.info$applyExpressionList)
  cat("\tapplyExpressionList duration:", (proc.time()-pt2)[3], "s\n")

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
  # TODO: set correct time values
  dimList[[(valueList[[1]] |> dimnames() |> names())[3]]] <- seq_len((valueList[[1]] |> dim())[3])

  pt3 <- proc.time()
  saveNetCdf(
    outFilePath,
    dimList = dimList,
    valueList = valueList)
  cat("\tSaving to", outFilePath, "duration:", (proc.time()-pt3)[3], "s\n")
  cat("All done after", (proc.time()-pt)[3], "s\n")
}


getDataAllTimeAll <- function(lotIdx, lineCount) {
  res <- lapply(names(.info$data), getDataTimeAll, lotIdx = lotIdx, lineCount = lineCount)
  names(res) <- names(.info$data)
  return(res)
}


rollEval <- function(expr, dataList, size = 3) {
  # TODO
  n <- dataList[[1]] |> length()
  res <- sapply((1+size):(n-size), \(i) {
    d <- lapply(dataList, \(x) x[(i-size):(i+size)])
    rlang::eval_tidy(expr, data = d)
  })
  return(c(rep(NA_real_, size), res, rep(NA_real_, size)))
}
