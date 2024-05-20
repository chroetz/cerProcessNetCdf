#' @export
concatAfterAggAllTime <- function(
  dirPath,
  pattern,
  outFilePath,
  deflate = 9
) {

  filePaths <- list.files(
    dirPath,
    pattern = pattern,
    full.names = TRUE)

  meta <- getNetCdfLatLonMeta(filePaths, latDecreasing = TRUE) # TODO: assumes lat order decreasing
  varNames <- meta[[1]]$varName
  lons <- meta[[1]]$lon
  lats <- lapply(meta, \(m) m$lat) |> unlist()
  stopifnot(order(lats) == seq_along(lats))

  dimList <- list(lon = lons, lat = lats)
  outNc <- initNetCdf(outFilePath, dimList, varNames, deflate, close=FALSE)

  for (info in meta) {
    cat("Processing lat range", info$lat[1], "to", info$lat[length(info$lat)], "...")
    pt <- proc.time()
    nc <- openNc(info$filePath)
    data <- read.nc(nc)
    close.nc(nc)
    latIdxs <- which(lats %in% info$lat)
    start <- c(NA, min(latIdxs))
    count <- c(NA, length(latIdxs))
    for (nm in varNames) {
      var.put.nc(outNc, nm, data[[nm]], start = start, count = count)
    }
    cat("done after", (proc.time()-pt)[3], "s\n")
  }
  close.nc(outNc)

  return(invisible(NULL))
}

getNetCdfLatLonMeta <- function(filePaths, latDecreasing) {

  meta <- lapply(filePaths, \(filePath) {
    nc <- openNc(filePath)
    res <- list(
      filePath = filePath,
      varName = ncGetNonDimVariableNames(nc),
      lon = var.get.nc(nc, "lon"),
      lat = var.get.nc(nc, "lat"))
    close.nc(nc)
    return(res)
  })

  m1 <- meta[[1]]
  for (m in meta) stopifnot(all(m$lon == m1$lon))
  latFirsts <- sapply(meta, \(m) m$lat[1])
  latLasts <- sapply(meta, \(m) m$lat[length(m$lat)])
  latMins <- sapply(meta, \(m) min(m$lat))
  latMaxs <- sapply(meta, \(m) max(m$lat))

  if (latDecreasing) {
    stopifnot(all(latFirsts == latMaxs))
    stopifnot(all(latLasts == latMins))
  } else {
    stopifnot(all(latFirsts == latMins))
    stopifnot(all(latLasts == latMaxs))
  }

  order <- order(latFirsts, decreasing = latDecreasing)
  meta <- meta[order]

  latFirstsNew <- sapply(meta, \(m) m$lat[1])
  stopifnot(order(latFirstsNew, decreasing = latDecreasing) == seq_along(latFirstsNew))

  return(meta)
}

