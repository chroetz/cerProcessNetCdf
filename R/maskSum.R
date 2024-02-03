#' @export
sumMask <- function(
  maskFilePath,
  outFilePath,
  regionRegex = NULL
) {

  # put function arguments into environment .info
  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(
    argNames,
    \(nm) assign(nm, env[[nm]], .info))

  gridFormat <- getNativeGridFormatFromFile(maskFilePath)
  initializeGrid(gridFormat)

  .info$targetFormat <- gridFormat
  openAndCheckMaskFile(maskFilePath)

  ptOuter <- proc.time()

  regionNames <- .info$maskList$regionNames
  if (hasValue(regionRegex)) regionNames <- str_subset(regionNames, regionRegex)
  cat(length(regionNames), "regions to process.\n")

  maskScalingValues <- NULL
  for (regionName in regionNames) {
    pt <- proc.time()
    cat("Process region '", regionName, "' ... ", sep="")
    maskValues <- getMaskValues(regionName, .info$maskList)
    if (is.null(maskScalingValues)) {
      maskScalingValues <- maskValues
    } else {
      maskScalingValues <- maskScalingValues + maskValues
    }
    cat("Done in ", (proc.time()-pt)[3], "s.\n")
  }

  cat("Close mask NC-File ... ")
  close.nc(.info$maskList$nc)
  cat("Done.\n")

  cat("get maskScalingValues duration:", (proc.time()-ptOuter)[3], "s\n")

  cat("Save maskScalingValues ... ")
  saveNetCdf(
    .info$outFilePath,
    dimList = list(lon = .info$grid$lonValues, lat = .info$grid$latValues),
    valueList = list(maskSum = maskScalingValues))
  cat("Done.\n")
}
