#' @export
setupSumAggregation <- function(
  targetFormat,
  maskFilePath,
  maskSumFilePath,
  boundingBoxFilePath,
  variableDataDescriptor,
  outFilePath
) {
  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(
    argNames,
    \(nm) assign(nm, env[[nm]], .info))

  initializeGrid(targetFormat)

  openAndCheckMaskFile(maskFilePath)
  readBoundingBoxes(boundingBoxFilePath)
  readMaskSum(maskSumFilePath)
  loadData("variable", variableDataDescriptor)

  return(invisible())
}


#' @export
runSumAggregation <- function(
    nBatches = 1,
    batchIndex = 1,
    yearsFilter = NULL,
    regionIndices = NULL
) {

  years <- getDataYears("variable")
  if (hasValue(yearsFilter)) years <- intersect(years, yearsFilter)
  cat(length(years), "years to process in total.\n")

  regionNames <- .info$maskList$regionNames
  if (hasValue(regionFilter)) regionNames <- intersect(regionNames, regionFilter)
  cat(length(regionNames), "regions to process.\n")

  cat("Split years into ", nBatches, "batches.\n")
  batches <- setupBatches(years, nBatches)
  batch <- batches[[batchIndex]]
  cat("Process batch", batchIndex, "with", length(batch), "years\n")

  cat("Start main loop.\n")
  for (year in batch) {
    cat("Starting year", year, "\n")
    ptYear <- proc.time()
    variableValuesAll <- getData("variable", year, setNaToZero = TRUE)
    scaledVariableValuesAll <- variableValuesAll / .info$maskSum$maskSum
    values <- vapply(
      regionNames,
      \(regionName) {
        pt <- proc.time()
        cat("\tRegion ", regionName, "... ")
        scaledVariableValuesRegion <- subsetRegion(
          scaledVariableValuesAll, regionName, .info$boundingBoxes, .info$grid)
        maskValues <- getMaskValues(regionName, .info$maskList, .info$boundingBoxes)
        value <- sum(maskValues * scaledVariableValuesRegion, na.rm = TRUE)
        cat("done after", (proc.time()-pt)[3], "s\n")
        return(value)
      },
      double(1))
    result <- tibble(
      year = year,
      region = regionNames,
      value = values)
    cat("Write year", year, "values to file", .info$outFilePath, "\n")
    readr::write_csv(result, .info$outFilePath, append = TRUE)
    cat("Year", year, "done after", (proc.time()-ptYear)[3], "s\n")
  }
  cat("End main loop.\n")

  cat("Close mask NC-File ... ")
  close.nc(.info$maskList$nc)
  cat("Done.\n")
}
