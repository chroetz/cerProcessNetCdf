hasValue <- function(x) {
  if (missing(x)) return(FALSE)
  return(length(x) > 0)
}


hasValueString <- function(x) {
  if (missing(x)) return(FALSE)
  if (length(x) == 0) return(FALSE)
  if (!is.character(x)) return(FALSE)
  if (is.na(x[1])) return(FALSE)
  if (nchar(x[1]) == 0) return(FALSE)
  return(TRUE)
}


removeFileNameEnding <- function(x) {
  return(gsub("\\.[^.]*$", "", x))
}
