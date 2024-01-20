#' @export
aggregateNaryMasked <- function(
  targetFormat = getNativeGridFormatFromFile(maskFilePath),
  maskFilePath,
  maskSumFilePath = NULL,
  boundingBoxFilePath = NULL,
  variableDataDescriptorList,
  aggregateExpression,
  outFilePath,
  outMetaFilePath = str_replace(outFilePath, "\\.csv$", "_meta.csv"),
  yearsFilter = NULL,
  regionRegex = NULL,
  nBatches = 1,
  batchIndex = 1
) {

  clearInfo()
  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(
    argNames,
    \(nm) assign(nm, env[[nm]], .info))

  initializeGrid(targetFormat)

  openAndCheckMaskFile(maskFilePath)
  if (hasValueString(boundingBoxFilePath)) readBoundingBoxes(boundingBoxFilePath)
  if (hasValueString(maskSumFilePath)) readMaskSum(maskSumFilePath)
  for (variableDataDescriptor in variableDataDescriptorList$list) {
    loadData(variableDataDescriptor)
  }

  regionNames <- .info$maskList$regionNames
  if (hasValue(regionRegex)) regionNames <- str_subset(regionNames, regionRegex)
  cat(length(regionNames), "regions to process.\n")

  labelsAndYears <- getDataLabelsAndYearsAll(yearsFilter)
  batch <- EbmUtility::splitAndGetOneBatch(
    "label-year-combinations",
    seq_len(nrow(labelsAndYears)),
    nBatches,
    batchIndex)

  if (file.exists(.info$outFilePath)) {
    cat(.info$outFilePath, "already exists. Deleting.\n")
    file.remove(.info$outFilePath)
  }
  if (file.exists(.info$outMetaFilePath)) {
    cat(.info$outMetaFilePath, "already exists. Deleting.\n")
    file.remove(.info$outMetaFilePath)
  }

  cat("Start main loop.\n")
  for (idx in batch) {
    labelsAndYear <- labelsAndYears[idx, ] |> as.list()
    labels <- labelsAndYear[-which(names(labelsAndYear) == "year")]
    processYearNaryAggregation(labels, labelsAndYear$year, regionNames)
  }
  cat("End main loop.\n")

  cat("Close mask NC-File ... ")
  close.nc(.info$maskList$nc)
  cat("Done.\n")
}


processYearNaryAggregation <- function(labels, year, regionNames) {
  cat("Processing labels", paste(labels, collapse=", "), "and year", year, "\n")
  ptYear <- proc.time()
  values <- vapply(
    regionNames,
    \(regionName) {
      pt <- proc.time()
      dataEnv <- rlang::env()
      cat("\tRegion ", regionName, "... ")
      bbInfo <- getSingleBoundingBox(.info$boundingBoxes, regionName)
      assignDataAll(year, labels, dataEnv, bbInfo)
      maskValues <- getMaskValues(regionName, .info$maskList, bbInfo)
      if (hasValue(.info$maskSum)) {
        maskSumValues <- subsetBox(.info$maskSum$maskSum, bbInfo) # TODO respect target format
        dataEnv$mask <- maskValues / maskSumValues
      } else {
        dataEnv$mask <- maskValues
      }
      value <- rlang::eval_tidy(
        .info$aggregateExpression,
        rlang::new_data_mask(dataEnv))
      cat("done after", (proc.time()-pt)[3], "s\n")
      return(value)
    },
    double(1))
  result <- tibble(
    label = label,
    year = year,
    region = regionNames,
    value = values)
  cat("Write values to file", .info$outFilePath, "\n")
  readr::write_csv(
    result,
    .info$outFilePath,
    append = TRUE,
    col_names = !file.exists(.info$outFilePath))
  cat("Year", year, "done after", (proc.time()-ptYear)[3], "s\n")
}
