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




ncLoadTimeDimension <- function(nc) {
  timeDimName <- setdiff(dimNames, c("lon", "lat"))
  stopifnot(length(timeDimName) == 1)
  timeValues <- var.get.nc(nc, timeDimName) |> as.vector()
  timeVarInfo <- var.inq.nc(nc, timeDimName)

  attNames <- sapply(
    seq_len(timeVarInfo$natts)-1,
    \(i) att.inq.nc(nc, timeDimName, i)$name)
  if ("units" %in% attNames) {
    timeUnitDescription <- att.get.nc(nc, timeDimName, "units")
    pattern <- "^days since ([\\d-]+)( \\d{2}:\\d{2}:(\\d{2})?)?"
    stopifnot(str_detect(timeUnitDescription, pattern))
    startDayText <- str_match(timeUnitDescription, pattern)[,2]
    startDate <- as.Date(startDayText)
    timeDates <- startDate + lubridate::days(timeValues)
    years <- lubridate::year(timeDates)
    formattedStartDate <- format(startDate, "%B %d, %Y")
    cat(
      "Assume that time values are days since year", startDate, "(", formattedStartDate, ").\n")
  } else {
    years <- timeValues
    cat("Assume that time values are years.\n")
  }

  if (max(abs(years - round(years))) >= sqrt(.Machine$double.eps)) {
    cat("PROBLEM in loadDataSingleFile with dataDescriptor\n")
    print(dataDescriptor)
    stop(
      "Could not correctly transform time dimension to years. Got:\n",
      paste0(years, collapse = ", "),
      "\nThese are not integer years.")
  }

  return(
    list(
      name = timeDimName,
      years = years,
      timeValues = timeValues))
}
