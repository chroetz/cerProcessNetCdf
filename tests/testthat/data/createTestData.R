library(RNetCDF)

createCountingFile <- function(n, years, fileName, dimOrder = 1:3, lonIncreasing = TRUE, latIncreasing = TRUE) {
  outPath <- testthat::test_path("data", fileName)
  nc <- create.nc(outPath, format = "netcdf4", share = FALSE)
  degStep <- 180/n
  dim.def.nc(nc, "lat", dimlength = n)
  dim.def.nc(nc, "year", dimlength = length(years))
  dim.def.nc(nc, "lon", dimlength = 2*n)
  var.def.nc(nc, "lon", "NC_DOUBLE", "lon")
  var.def.nc(nc, "lat", "NC_DOUBLE", "lat")
  var.def.nc(nc, "year", "NC_INT", "year")
  var.def.nc(nc, "value", "NC_DOUBLE", c("lon", "lat", "year")[dimOrder], deflate = 9)
  lonValues <- seq(-180, 180, by = degStep)[-1] - degStep/2
  if (!lonIncreasing) lonValues <- lonValues |> rev()
  latValues <- seq(-90, 90, by = degStep)[-1] - degStep/2
  if (!latIncreasing) latValues <- latValues |> rev()
  var.put.nc(nc, "lon", lonValues)
  var.put.nc(nc, "lat", latValues)
  var.put.nc(nc, "year", years)
  dataArray <- array(seq_len(2*n*n*length(years)), dim = c(2*n, n, length(years)))
  if (!lonIncreasing) dataArray <- dataArray[(2*n):1, , ]
  if (!latIncreasing) dataArray <- dataArray[, n:1, ]
  dataArray <- aperm(dataArray, dimOrder)
  var.put.nc(nc, "value", dataArray)
  close.nc(nc)
}

years <- 2000:2002
for (n in 2^(1:4)) {
  createCountingFile(n, years, sprintf("test_counting_%d_%d_%d.nc", n, min(years), max(years)))
}

n <- 4
for (dimOrder in list(
    c(1, 2, 3),
    c(1, 3, 2),
    c(2, 1, 3),
    c(2, 3, 1),
    c(3, 1, 2),
    c(3, 2, 1))) {
      for (lonIncreasing in c(TRUE, FALSE)) {
        for (latIncreasing in c(TRUE, FALSE)) {
          createCountingFile(
            n,
            years,
            sprintf(
              "test_counting_%d_%d_%d_%d_%d_%d_%s_%s.nc",
              n, min(years), max(years), dimOrder[1], dimOrder[2], dimOrder[3],
              if (lonIncreasing) "lonInc" else "lonDec",
              if (latIncreasing) "latInc" else "latDec"),
            dimOrder = dimOrder,
            lonIncreasing = lonIncreasing,
            latIncreasing = latIncreasing)
        }
      }
}
