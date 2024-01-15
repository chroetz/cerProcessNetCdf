hasValue <- function(x) {
  if (missing(x)) return(FALSE)
  return(length(x) > 0)
}


removeFileNameEnding <- function(x) {
  return(gsub("\\.[^.]*$", "", x))
}
