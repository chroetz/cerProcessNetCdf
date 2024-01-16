testLoadAndGet <- function(
    n, year, lonFirst,
    lonIncreasing, latIncreasing,
    inLonIncreasing, inLatIncreasing,
    dimOrder = 1:3
) {
  targetFormat <- ConfigOpts::makeOpts(
    "GridFormat",
    nLon = 2*n,
    nLat = n,
    lonIncreasing = lonIncreasing,
    latIncreasing = latIncreasing,
    lonFirst = lonFirst)
  initializeGrid(targetFormat)
  loadData("testData", ConfigOpts::makeOpts(c("SingleFile", "DataDescriptor"),
    filePath = test_path(
      "data",
      sprintf(
        "test_counting_%d_2000_2002_%d_%d_%d_%s_%s.nc",
        n, dimOrder[1], dimOrder[2], dimOrder[3],
        if (inLonIncreasing) "lonInc" else "lonDec",
        if (inLatIncreasing) "latInc" else "latDec"))))
  res <- getData("testData", year)
  expectedRes <- matrix(seq_len(2*n*n), nrow = 2*n)
  dimnames(expectedRes) <- list(lon = NULL, lat = NULL)
  if (!lonIncreasing) expectedRes <- expectedRes[(2*n):1, ]
  if (!latIncreasing) expectedRes <- expectedRes[, n:1]
  if (!lonFirst) expectedRes <- t(expectedRes)
  expect_equal(res, expectedRes)
}

test_that("transform to traget", {
  for (dimOrder in list(
    c(1, 2, 3),
    c(1, 3, 2),
    c(2, 1, 3),
    c(2, 3, 1),
    c(3, 1, 2),
    c(3, 2, 1))) {
    for (lonFirst in c(TRUE, FALSE)) {
      for (lonIncreasing in c(TRUE, FALSE)) {
        for (latIncreasing in c(TRUE, FALSE)) {
          for (inLonIncreasing in c(TRUE, FALSE)) {
            for (inLatIncreasing in c(TRUE, FALSE)) {
              testLoadAndGet(
                n = 4,
                year = 2000,
                lonFirst = lonFirst,
                lonIncreasing = lonIncreasing,
                latIncreasing = latIncreasing,
                inLonIncreasing = inLonIncreasing,
                inLatIncreasing = inLatIncreasing,
                dimOrder = dimOrder)
            }
          }
        }
      }
    }
  }
})
