helmert <- function(d){
  stopifnot(d >= 2L)
  mc <- contr.helmert(d)
  k <- seq_len(d - 1L)
  m <- sweep(mc, 2L, sqrt(k * (k + 1)), "/")
  t(m)
}