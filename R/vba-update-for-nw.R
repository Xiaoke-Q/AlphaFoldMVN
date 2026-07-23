.afmvn_expected_stats <- function(M,B,r) {
  N <- dim(M)[1]
  G <- dim(M)[2]
  p <- dim(M)[3]

  M_sum <- numeric(p)
  B_sum <- matrix(0, nrow = p, ncol = p)

  for (i in 1:N) {
    for (g in 1:2) {
      M_sum <- M_sum + r[i, g]*M[i, g, ]
      B_sum <- B_sum + r[i, g]*B[i, g, , ]
    }
  }

  list(
    M_sum = M_sum,
    B_sum = B_sum
  )
}

.afmvn_update_nw <- function(M, B, r,
                             m_0, kappa_0, Psi_0, nu_0) {

  stats <- .afmvn_expected_stats(M = M, B = B, r = r)

  N <- dim(M)[1]
  M_q <- stats$M_sum
  B_q <- stats$B_sum

  kappa_q <- kappa_0 + N
  nu_q <- nu_0 + N
  m_q <- (kappa_0*m_0 + M_q) / kappa_q

  Psi_q_inv <- solve(Psi_0) + B_q + kappa_0 * tcrossprod(m_0) - kappa_q * tcrossprod(m_q)
  Psi_q <- solve(Psi_q_inv)
  Psi_q <- 0.5 * (Psi_q + t(Psi_q))

  list(
    kappa_q = kappa_q,
    nu_q = nu_q,
    m_q = m_q,
    Psi_q = Psi_q,
    Psi_q_inv = Psi_q_inv
  )
}