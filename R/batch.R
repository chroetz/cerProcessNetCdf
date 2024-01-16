setupBatches <- function(vec, nBatches) {
  n <- length(vec)
  batchSizeSmall <- floor(n / nBatches)
  batchSizeLarge <- batchSizeSmall + 1
  nBatchesLarge <- n - batchSizeSmall * nBatches
  nBatchesSmall <- nBatches - nBatchesLarge
  batches <- split(
    vec,
    c(rep(seq_len(nBatchesLarge), each=batchSizeLarge),
      rep(nBatchesLarge+seq_len(nBatchesSmall), each=batchSizeSmall)))
  return(batches)
}
