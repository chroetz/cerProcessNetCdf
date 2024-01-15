#' @export
createDataDescriptorYearlyFiles <- function(dir, pattern) {
  dataDescriptor <- lst(
    type = "YerlyFiles",
    dir, pattern)
  class(dataDescriptor) <- c("DataDescriptor", "YearlyFiles")
  return(dataDescriptor)
}

#' @export
createDataDescriptorSingleFile <- function(filePath) {
  dataDescriptor <- lst(
    type = "SingleFile",
    filePath)
  class(dataDescriptor) <- c("DataDescriptor", "SingleFile")
  return(dataDescriptor)
}
