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
runSumAggregation <- function(yearsFilter = NULL, regionIndices = NULL) {

  years <- getDataYears("variable")
  if (hasValue(yearsFilter)) years <- intersect(years, yearsFilter)
  cat(length(years), "years to process.\n")

  regionNames <- .info$maskList$regionNames
  if (hasValue(regionIndices)) regionNames <- regionNames[regionIndices]
  cat(length(regionNames), "regions to process.\n")

  cat("Start main loop.\n")
  for (year in years) {
    cat("Year:", year, "\n")
    variableValuesAll <- getData("variable", year, setNaToZero = TRUE)
    scaledVariableValuesAll <- variableValuesAll / .info$maskSum$maskSum
    values <- vapply(
      regionNames,
      \(regionName) {
        pt <- proc.time()
        cat("\tRegion:", regionName, "\n")
        scaledVariableValuesRegion <- subsetRegion(
          scaledVariableValuesAll, regionName, .info$boundingBoxes, .info$grid)
        maskValues <- getMaskValues(regionName, .info$maskList, .info$boundingBoxes)
        value <- sum(maskValues * scaledVariableValuesRegion, na.rm = TRUE)
        cat("\tprocessRegionYear duration:", (proc.time()-pt)[3], "s\n")
        return(value)
      },
      double(1))
    result <- tibble(
      year = year,
      region = regionNames,
      value = values)
    readr::write_csv(result, .info$outFilePath, append = TRUE)
  }
  cat("End main loop.\n")

  cat("Close mask NC-File ... ")
  close.nc(.info$maskList$nc)
  cat("Done.\n")
}
