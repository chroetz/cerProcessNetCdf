#' @export
setupMaskSummation <- function(
  maskFilePath,
  outFilePath
) {

  # put function arguments into environment .info
  argNames <- rlang::fn_fmls_names()
  env <- rlang::current_env()
  lapply(
    argNames,
    \(nm) assign(nm, env[[nm]], .info))

  gridFormat <- getNativeGridFormatFromFile(maskFilePath)
  initializeGrid(gridFormat)

  openAndCheckMaskFile(maskFilePath)

  return(invisible())
}


#' @export
runMaskSummation <- function() {

  pt <- proc.time()

  regionNames <- .info$maskList$regionNames
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

  cat("getMaskScaling duration:", (proc.time()-pt)[3], "s\n")

  cat("Save mask scaling values ... ")
  saveNetCdf(
    .info$outFilePath,
    list(lon = .info$grid$lonValues, lat = .info$grid$latValues),
    maskScalingValues)
  cat("Done.\n")
}
