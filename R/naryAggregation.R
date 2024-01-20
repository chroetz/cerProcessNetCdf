#' @export
aggregateNaryMasked <- function(
  targetFormat = getNativeGridFormatFromFile(maskFilePath),
  maskFilePath,
  maskSumFilePath = NULL,
  boundingBoxFilePath = NULL,
  variableDataDescriptorList,
  aggregateExpression = rlang::expr(sum(mask * Reduce(`*`, lapply(ls(), get)), na.rm=TRUE)),
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
  if (file.exists(.info$outMetaFilePath)) {
    cat(.info$outMetaFilePath, "already exists. Deleting.\n")
    file.remove(.info$outMetaFilePath)
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
      dataEnv <- rlang::env()
      cat("\tRegion ", regionName, "... ")
      bbInfo <- getSingleBoundingBox(.info$boundingBoxes, regionName)
      assignDataAll(year, dataEnv, bbInfo)
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
    year = year,
    region = regionNames,
    value = values)
  cat("Write year", year, "values to file", .info$outFilePath, "\n")
  readr::write_csv(
    result,
    .info$outFilePath,
    append = TRUE,
    col_names = !file.exists(.info$outFilePath))
  writeInfo(
    year,
    .info$outMetaFilePath,
    append = TRUE,
    col_names = !file.exists(.info$outMetaFilePath))
  cat("Year", year, "done after", (proc.time()-ptYear)[3], "s\n")
}


writeInfo <- function(year, outFilePath, ...) {
  meta <- lapply(names(.info$data), getInfo, year = year) |> bind_rows()
  readr::write_csv(meta, outFilePath, ...)
}
