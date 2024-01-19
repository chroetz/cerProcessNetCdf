#' @export
aggregateNaryMasked <- function(
  targetFormat = getNativeGridFormatFromFile(maskFilePath),
  maskFilePath,
  maskSumFilePath = NULL,
  boundingBoxFilePath = NULL,
  variableDataDescriptorList,
  aggregateFunction = \(mask, variables) sum(mask * Reduce(`*`, variables), na.rm=TRUE),
  outFilePath,
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
  for (variableDataDescriptor in variableDataDescriptorList) {
    loadData(variableDataDescriptor)
  }

  years <- getDataYearsAll()
  if (hasValue(yearsFilter)) years <- intersect(years, yearsFilter)
  cat(
    length(years),
    "years to process in total, from",
    min(years),
    "to",
    max(years),
    "\n")

  regionNames <- .info$maskList$regionNames
  if (hasValue(regionRegex)) regionNames <- str_subset(regionNames, regionRegex)
  cat(length(regionNames), "regions to process.\n")

  batch <- EbmUtility::splitAndGetOneBatch("years", years, nBatches, batchIndex)

  if (file.exists(.info$outFilePath)) {
    cat(.info$outFilePath, "already exists. Deleting.\n")
    file.remove(.info$outFilePath)
  }

  cat("Start main loop.\n")
  for (year in batch) {
    processYearNaryAggregation(year, regionNames)
  }
  cat("End main loop.\n")

  cat("Close mask NC-File ... ")
  close.nc(.info$maskList$nc)
  cat("Done.\n")
}


processYearNaryAggregation <- function(year, regionNames) {
  cat("Processing year", year, "\n")
  ptYear <- proc.time()
  values <- vapply(
    regionNames,
    \(regionName) {
      pt <- proc.time()
      cat("\tRegion ", regionName, "... ")
      bbInfo <- getSingleBoundingBox(.info$boundingBoxes, regionName)
      variablesValuesList <- getDataAll(year, bbInfo)
      maskValues <- getMaskValues(regionName, .info$maskList, bbInfo)
      if (hasValue(.info$maskSum)) {
        maskSumValues <- subsetBox(.info$maskSum$maskSum, bbInfo) # TODO respect target format
        value <- .info$aggregateFunction(maskValues / maskSumValues, variablesValuesList)
      } else {
        value <- .info$aggregateFunction(maskValues, variablesValuesList)
      }
      cat("done after", (proc.time()-pt)[3], "s\n")
      return(value)
    },
    double(1))
  result <- tibble(
    year = year,
    region = regionNames,
    value = values)
  cat("Write year", year, "values to file", .info$outFilePath, "\n")
  readr::write_csv(
    result,
    .info$outFilePath,
    append = TRUE,
    col_names = !file.exists(.info$outFilePath))
  cat("Year", year, "done after", (proc.time()-ptYear)[3], "s\n")
}
