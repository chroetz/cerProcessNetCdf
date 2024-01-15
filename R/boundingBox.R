#' @export
saveBoundingBoxes <- function(
    boundingBoxes,
    outFilePath,
    maskFilePath,
    regionVariableName = "regionName"
  ) {

  maskNc <- open.nc(maskFilePath)
  on.exit(close.nc(maskNc))

  outNc <- create.nc(outFilePath, format = "netcdf4", share = FALSE)
  on.exit(close.nc(outNc))

  # copy dimensions from maskNc to outNc
  for (i in seq_len(file.inq.nc(maskNc)$ndims)) {
    dimInfo <- dim.inq.nc(maskNc, i - 1)
    varInfo <- var.inq.nc(maskNc, dimInfo$name)
    dim.def.nc(outNc, dimInfo$name, dimInfo$length)
    var.def.nc(outNc, dimInfo$name, varInfo$type, dimInfo$name)
    var.put.nc(outNc, dimInfo$name, var.get.nc(maskNc, dimInfo$name))
    for (j in seq_len(varInfo$natts)) {
      attInfo <- att.inq.nc(maskNc, dimInfo$name, j - 1)
      att.put.nc(
        outNc,
        dimInfo$name,
        attInfo$name,
        attInfo$type,
        att.get.nc(maskNc, dimInfo$name, attInfo$name))
    }
  }

  # create regions dimension
  regions <- colnames(boundingBoxes)
  dim.def.nc(outNc, regionVariableName, length(regions))
  var.def.nc(outNc, regionVariableName, "NC_STRING", regionVariableName)
  var.put.nc(outNc, regionVariableName, regions)

  # put bounding boxes values into outNc
  for (i in seq_len(nrow(boundingBoxes))) {
    varName <- rownames(boundingBoxes)[i]
    var.def.nc(outNc, varName, "NC_DOUBLE", regionVariableName)
    var.put.nc(outNc, varName, boundingBoxes[i, ])
  }
}

#' @export
getBoundingBoxesFromMask <- function(path) {
  nc <- open.nc(path)
  on.exit(close.nc(nc))
  fileInfo <- file.inq.nc(nc)
  stopifnot(fileInfo$ndims == 2)
  ndims <- 2
  dimNames <- c(dim.inq.nc(nc, 0)$name, dim.inq.nc(nc, 1)$name)
  boundingBoxes <- sapply(
    seq_len(fileInfo$nvars - ndims),
    \(i) {
      cat("Processing region ", i, "... ")
      pt <- proc.time()
      mask <- var.get.nc(nc, ndims + i - 1)
      bbox <- getBoundingBox(mask)
      cat("done in ", (proc.time() - pt)[3], "s.\n")
      return(bbox)
    })
  rownames(boundingBoxes) <- c(
    paste0("min_", dimNames[1]),
    paste0("max_", dimNames[1]),
    paste0("min_", dimNames[2]),
    paste0("max_", dimNames[2]))
  colnames(boundingBoxes) <- sapply(
    seq_len(fileInfo$nvars-2),
    \(i) var.inq.nc(nc, 1 + i)$name)
  return(boundingBoxes)
}

getBoundingBox <- function(x) {
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(x))
  rowSumsX <- x |> abs() |> rowSums()
  rowBoundingIdxs <- range(which(rowSumsX > 0))
  colSumsX <- x |> abs() |> colSums()
  colBoundingIdxs <- range(which(colSumsX > 0))
  return(c(
    rowMin = rowBoundingIdxs[1],
    rowMax = rowBoundingIdxs[2],
    colMin = colBoundingIdxs[1],
    colMax = colBoundingIdxs[2]))
}
