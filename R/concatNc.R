#' @export
runConcatNetCdf <- function(
  outFilePath,
  inDirPath,
  inFilePattern
) {

  inFileNames <- dir(inDirPath, pattern = inFilePattern)
  cat("Found", length(inFileNames), "files.\n")

  # Assume all files have the same dims (in the same order).
  firstInNc <- openNc(file.path(inDirPath, inFileNames[1]))
  fileInfo <- file.inq.nc(firstInNc)
  cat("Create output file.\n")
  outNc <- create.nc(outFilePath, format = "netcdf4", share = FALSE, prefill = FALSE)

  # copy dimensions from firstInNc to outNc
  cat("Copy dimensions from file", inFileNames[1], "to output file.\n")
  for (i in seq_len(fileInfo$ndims)) {
    dimInfo <- dim.inq.nc(firstInNc, i - 1)
    varInfo <- var.inq.nc(firstInNc, dimInfo$name)
    dim.def.nc(outNc, dimInfo$name, dimInfo$length)
    var.def.nc(outNc, dimInfo$name, varInfo$type, dimInfo$name)
    var.put.nc(outNc, dimInfo$name, var.get.nc(firstInNc, dimInfo$name))
    for (j in seq_len(varInfo$natts)) {
      attInfo <- att.inq.nc(firstInNc, dimInfo$name, j - 1)
      att.put.nc(
        outNc,
        dimInfo$name,
        attInfo$name,
        attInfo$type,
        att.get.nc(firstInNc, dimInfo$name, attInfo$name))
    }
  }

  close.nc(firstInNc)

  cat("Copy variables from all files to output file.\n")
  for (fileName in inFileNames) {

    pt <- proc.time()
    cat("Processing", fileName, "\n")
    inNc <- openNc(file.path(inDirPath, fileName))
    varNames <- ncGetNonDimVariableNames(inNc)

    # copy variables from inNc to outNc (requires dims to have the same order as in firstNc)
    for (varName in varNames) {
      cat(varName, ",", sep="")
      varInfo <- var.inq.nc(inNc, varName)
      var.def.nc(outNc, varName, varInfo$type, varInfo$dimids, deflate=9)
      var.put.nc(outNc, varName, var.get.nc(inNc, varName))
      for (j in seq_len(varInfo$natts)) {
        attInfo <- att.inq.nc(inNc, varName, j - 1)
        att.put.nc(
          outNc,
          varName,
          attInfo$name,
          attInfo$type,
          att.get.nc(inNc, varName, attInfo$name))
      }
    }
    cat("\n")
    close.nc(inNc)
    cat("Processing", fileName, "took", (proc.time() - pt)[3], "s.\n")
  }

  cat("Close output file.\n")
  close.nc(outNc)
  cat("Done.\n")
}
