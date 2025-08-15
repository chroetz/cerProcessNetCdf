getFileTimes <- function(ncFilePath, timeDimName) {
  nc <- openNc(ncFilePath)
  times <- getNcTimes(nc, timeDimName)
  close.nc(nc)
  return(times)
}

getNcTimes <- function(nc, timeDimName) {
  values <- varGetNc(nc, timeDimName)
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
    seconds = origin + lubridate::seconds(values),
    minutes = origin + lubridate::minutes(values),
    hours = origin + lubridate::hours(values),
    days = origin + lubridate::days(values),
    weeks = origin + lubridate::weeks(values),
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




ncLoadTimeDimension <- function(nc, timeDimName = NULL) {

  # TODO: Assumes calender proleptic_gregorian. Warn if not.

  if (is.null(timeDimName)) {
    dimNames <- ncGetDimensionNames(nc)
    timeDimName <- setdiff(dimNames, c("lon", "lat"))
  }

  stopifnot(length(timeDimName) == 1)
  timeValues <- varGetNc(nc, timeDimName) |> as.vector()
  timeVarInfo <- var.inq.nc(nc, timeDimName)

  attNames <- sapply(
    seq_len(timeVarInfo$natts)-1,
    \(i) att.inq.nc(nc, timeDimName, i)$name)
  if ("units" %in% attNames) {
    timeUnitDescription <- att.get.nc(nc, timeDimName, "units")
    patternDaysSince <- "^days since ([\\d-]+)( \\d{2}:\\d{2}:(\\d{2})?)?"
    patternYearsSince <- "^years since ([\\d]+)"
    if (str_detect(timeUnitDescription, patternDaysSince)) {
      startDayText <- str_match(timeUnitDescription, patternDaysSince)[,2]
      startDate <- as.Date(startDayText)
      timeDates <- startDate + lubridate::days(timeValues)
      years <- lubridate::year(timeDates)
      formattedStartDate <- format(startDate, "%B %d, %Y")
      cat("Assume that time values are days since year", startDate, "(", formattedStartDate, ").\n")
    } else if (str_detect(timeUnitDescription, patternYearsSince)) {
      startYearText <- str_match(timeUnitDescription, patternYearsSince)[,2]
      startYear <- as.integer(startYearText)
      stopifnot(is.finite(startYear))
      years <- startYear + timeValues
      cat("Assume that time values are years since year", startYear, ".\n")
    } else {
      stop("Unknown time format: ", timeUnitDescription)
    }
  } else {
    years <- timeValues
    cat("Assume that time values are years.\n")
  }

  if (max(abs(years - round(years))) <= 1/360) { # tolerance for year value extraction of about one day
    years <- round(years)
  } else {
    warning(
      "Could not perfectly transform time dimension to years. Got:\n",
      paste0(years, collapse = ", "),
      "\nThese are not integer years. Rounding them to integer. ",
      "Description of time unit is ", timeUnitDescription)
    years <- round(years)
  }

  return(
    list(
      name = timeDimName,
      years = years,
      timeValues = timeValues))
}
