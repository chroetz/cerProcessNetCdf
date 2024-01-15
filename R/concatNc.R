#' @export
runConcatNetCdf <- function(
  outFilePath,
  inFileDir,
  inFilePattern
) {

  inFileNames <- dir(inFileDir, pattern = inFilePattern)
  cat("Found", length(inFileNames), "files.\n")

  # Assume all files have the same dims (in the same order).
  firstInNc <- open.nc(file.path(inFileDir, inFileNames[1]))
  fileInfo <- file.inq.nc(firstInNc)
  cat("Create output file.\n")
  outNc <- create.nc(outFilePath, format = "netcdf4")
  on.exit(close.nc(outNc))

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

    cat("Processing", fileName, "\n")
    inNc <- open.nc(file.path(inFileDir, inFileNames[1]))
    varNames <- ncGetNonDimVariableNames(inNc)

    # copy variables from inNc to outNc (requires dims to have the same order as in firstNc)
    for (varName in varNames) {
      varInfo <- var.inq.nc(inNc, varName)
      var.def.nc(outNc, varName, varInfo$type, varInfo$dimids)
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

    close.nc(inNc)
  }

  cat("Done.\n")
}