getFileTimes <- function(ncFilePath, timeDimName) {
  nc <- open.nc(ncFilePath)
  times <- getNcTimes(nc, timeDimName)
  close.nc(nc)
  return(times)
}

getNcTimes <- function(nc, timeDimName) {
  values <- var.get.nc(nc, timeDimName)
  unitText <- att.get.nc(nc, timeDimName, "units")
  times <- numericToTime(values, unitText)
  return(times)
}

numericToTime <- function(values, unitText) {
  stopifnot(is.numeric(values))
  stopifnot(is.character(unitText))
  stopifnot(length(unitText) == 1)
  match <- str_match(unitText, "(.*) since (.*)")
  stopifnot(!any(is.na(match)))
  origin <- lubridate::as_date(match[3])
  times <- switch(match[2],
    days = origin + lubridate::days(values),
    months = origin + lubridate::months(values),
    years = origin + lubridate::years(values),
    stop("Unknown time unit: ", match[2])
  )
  return(times)
}


getTimeIndicesYear <- function(ncFilePath, timeDimName, year) {
  times <- getFileTimes(ncFilePath, timeDimName)
  indices <- which(lubridate::year(times) == year)
  return(indices)
}
