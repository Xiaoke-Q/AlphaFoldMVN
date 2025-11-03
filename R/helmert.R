#' Generate an orthonormal Helmertian matrix in the strict sense, without the first row.
#'
#' @param d The number of columns of the matrix.
#'
#' @returns A (d-1)-by-d Helmertian matrix.
#'
#' @export
#' @examples
#' helmert(d = 4)
helmert <- function(d){
  stopifnot(d >= 2L)
  mc <- contr.helmert(d)
  k <- seq_len(d - 1L)
  m <- sweep(mc, 2L, sqrt(k * (k + 1)), "/")
  t(m)
}