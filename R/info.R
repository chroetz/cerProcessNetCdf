#' @export
getNcInfo <- function(filePath) {
  nc <- open.nc(filePath)
  ncInfo <- file.inq.nc(nc)
  dimInfo <- lapply(seq_len(ncInfo$ndims)-1, dim.inq.nc, ncfile = nc)
  varInfo <- lapply(seq_len(ncInfo$nvars)-1, var.inq.nc, ncfile = nc)
  attInfo <- lapply(
    varInfo,
    \(vi) {
      lapply(
        seq_len(vi$natts)-1,
        att.inq.nc,
        ncfile = nc,
        variable = vi$id)
    })
  attValue <- lapply(
    varInfo,
    \(vi) {
      lapply(
        seq_len(vi$natts)-1,
        att.get.nc,
        ncfile = nc,
        variable = vi$id)
    })
  close.nc(nc)
  return(lst(ncInfo, dimInfo, varInfo, attInfo, attValue))
}
