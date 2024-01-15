#' Setup for runPopWeightAggregation().
#' @param targetFormat An object of class GridFormat as created by createGridFormat().
#' @export
setupPopWeightAggregation <- function(
  targetFormat,
  maskFilePath,
  maskSumFilePath,
  boundingBoxFilePath,
  popDataDescriptor = NULL,
  invarDir,
  invarFileNamePattern,
  invarDimensionName,
  invarValueVariableName,
  lonLatVarToDimOrder = 1:3, # TODO: infer from file
  outDir,
  outNcFilePattern,
  batchSize
) {
  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(
    argNames,
    \(nm) assign(nm, env[[nm]], .info)
  )

  initializeGrid(targetFormat)

  openAndCheckMaskFile(maskFilePath)
  readBoundingBoxes(boundingBoxFilePath)
  readMaskSum(maskSumFilePath)

  if (hasValue(popDataDescriptor)) {
    loadData("population", popDataDescriptor)
    .info$weightByPop <- TRUE
  } else {
    .info$weightByPop <- FALSE
  }

  invarFileNames <- list.files(
    invarDir,
    pattern = invarFileNamePattern)
  fileYears <- stringr::str_match(invarFileNames, invarFileNamePattern)[,2] |> as.integer()
  .info$invarFileMeta <- tibble(
    year = fileYears,
    name = invarFileNames,
    path = file.path(invarDir, invarFileNames))
  return(invisible())
}


#' @export
setupPopWeightAggregationStatistics <- function(...) {
  args <- list(...)
  .info$statisticNames <- names(args)
  .info$statisticFunctions <- args
  return(invisible())
}


#' @export
runPopWeightAggregation <- function(
    yearsFilter = NULL,
    invarNamesIdxFilter = NULL,
    regionIndices = NULL
) {

  years <- getYearsPop()
  if (hasValue(yearsFilter)) years <- intersect(years, yearsFilter)
  cat(length(years), "years to process.\n")

  regionNames <- .info$maskList$regionNames
  if (hasValue(regionIndices)) regionNames <- regionNames[regionIndices]
  cat(length(regionNames), "regions to process.\n")

  cat("Initializing output files...\n")
  initOutNc(years, regionNames, .info$statisticNames)
  cat("Done.\n")

  cat("Start main loop.\n")
  for (year in years) {
    cat("Year:", year, "\n")

    outNcFilePath <- getOutNcFilePath(year)
    outNc <- open.nc(outNcFilePath, write = TRUE, share = FALSE)

    if (.info$weightByPop) {
      popValuesAll <- getData("population", year, setNaToZero = TRUE)
    }
    invarNames <- getInvarNames(year)
    if (hasValue(invarNamesIdxFilter)) {
      invarNames <- invarNames[invarNamesIdxFilter]
      invarNames <- invarNames[!is.na(invarNames)]
    }
    if (length(invarNames) == 0) {
      cat("No invarNames to process in year", year, ". Skipping.\n")
      next
    }
    checkInvar(year, invarNames)
    fullyFilledRegionNames <- getFullyFilledRegionNames(.info, year, invarNames, outNc=outNc)
    if (length(fullyFilledRegionNames) > 0) {
      cat(
        "\tFound",
        length(fullyFilledRegionNames),
        "regions with data. Not re-calculating those.\n")
      regionNames <- setdiff(regionNames, fullyFilledRegionNames)
    } else {
      cat("\tNo filled regions found. Processing all.\n")
    }
    for (regionName in regionNames) {
      cat("\tRegion:", regionName, "\n")
      pt <- proc.time()
      maskValues <- getMaskValues(regionName, .info$maskList, .info$idxBoundingBoxes)
      if (.info$weightByPop) {
        popValuesRegion <- subsetRegion(popValuesAll, regionName, .info$idxBoundingBoxes, .info$grid)
        aggregationDistri <- calculateProductDistribution(popValuesRegion, maskValues)
      } else {
        aggregationDistri <- normalizeDistribution(maskValues)
      }
      cat("\tload mask duration:", (proc.time()-pt)[3], "s\n")
      pt <- proc.time()
      processRegionYear(
        regionName,
        year,
        invarNames,
        aggregationDistri,
        batchSize = .info$batchSize,
        outNc = outNc)
      cat("\tprocessRegionYear duration:", (proc.time()-pt)[3], "s\n")
    }

    close.nc(outNc)
  }
  cat("End main loop.\n")

  cat("Close mask NC-File ... ")
  close.nc(.info$maskList$nc)
  cat("Done.\n")
}


getYearsPop <- function() {
  if (.info$weightByPop) {
    return(
      intersect(
        .info$data$population$years,
        .info$invarFileMeta$year))
  } else {
    return(.info$invarFileMeta$year)
  }
}
