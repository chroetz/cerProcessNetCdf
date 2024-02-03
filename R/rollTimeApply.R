rollTimeApply <- function(
  targetFormat = NULL,
  variableDataDescriptorList,
  applyExpressionList,
  outFilePath,
  pastTimeSteps = 0,
  futureTimeSteps = 0,
  paddingValue = NA_real_,
  paddingMode = c("crop", "extend", "value"),
  yearsFilter = NULL,
  nBatches = 1,
  batchIndex = 1
) {

  paddingMode <- match.arg(paddingMode)

  for (variableDataDescriptor in variableDataDescriptorList$list) {
    loadData(variableDataDescriptor)
  }

  if (is.null(targetFormat)) targetFormat <- .info$data[[1]]$gridFormat

  clearInfo()
  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(argNames, \(nm) assign(nm, env[[nm]], .info))

  initializeGrid(targetFormat)

  if (file.exists(.info$outFilePath)) {
    cat(.info$outFilePath, "already exists. Deleting.\n")
    file.remove(.info$outFilePath)
  }

  # Split grid into nBatches areas of roughly equal size.
  # batch <- cerUtility::splitAndGetOneBatch(
  #   "label-year-combinations",
  #   seq_len(.info$grid$nLon),
  #   nBatches,
  #   batchIndex)

  cat("Start main loop.\n")
  for (idx in batch) { # TODO
    processRollTimeApply(idx)
  }
  cat("End main loop.\n")
}


processRollTimeApply <- function(patchIndex) {

  # TODO: process one area patch for all times and save it

  cat("Processing labels", paste(labels, collapse=", "), "and year", year, "\n")
  ptYear <- proc.time()

  dataEnv <- rlang::env()
  assignDataAllRollTime(
    year, labels, dataEnv,
    pastTimeSteps,
    futureTimeSteps,
    paddingValue,
    paddingMode)
  dataMask <- rlang::new_data_mask(dataEnv)
  # TODO: time dimension:
  # load full time series but only a local patch of the grid
  if (.info$grid$lonFirst) {
    proto <- matrix(NA_real_, nrow = .info$grid$nLon, ncol = .info$grid$nLat)
    dimList <- list(lon = .info$grid$lonValues, lat = .info$grid$latValues)
  } else {
    proto <- matrix(NA_real_, nrow = .info$grid$nLat, ncol = .info$grid$nLon)
    dimList <- list(lat = .info$grid$lonValues, lon = .info$grid$latValues)
  }
  valueList <- lapply(
    .info$applyExpressionList,
    rlang::eval_tidy,
    data = dataMask,
    proto)
  names(valueList) <- names(.info$applyExpressionList)

  saveNetCdf(
    .info$outFilePath,
    dimList = dimList,
    valueList = valueList)

  cat(
    "Labels", paste(labels, collapse=", "),
    "and year", year,
    "done after", (proc.time()-ptYear)[3], "s\n")
}
