.afmvn_nw_stats <- function(data_ctx, gh_ctx, 
                            m_eta, s2_eta, r,
                            eps_alpha = 1e-8) {
  p <- data_ctx$p

  q <- .eta_quadrature(m_eta = m_eta, s2_eta = s2_eta, gh_ctx = gh_ctx)

  M_sum <- numeric(p)
  B_sum <- matrix(0, nrow = p, ncol = p)

  for (ell in seq_along(q$eta)) {
    comp <- .afmvn_components_eta(data_ctx = data_ctx, eta = q$eta[ell], eps_alpha = eps_alpha)
    w <- q$weight[ell]

    z <- comp$z0
    rg <- r[, 1]

    M_sum <- M_sum + w * colSums(z * rg)
    B_sum <- B_sum + w * crossprod(z * sqrt(rg))

    if (!isTRUE(comp$near_zero)) {
      z <- comp$z1
      rg <- r[, 2]

      M_sum <- M_sum + w * colSums(z * rg)
      B_sum <- B_sum + w * crossprod(z * sqrt(rg))
    }
  }

  list(M_sum = M_sum, B_sum = B_sum)
}

.afmvn_update_nw <- function(M_sum, B_sum, N,
                             m_0, kappa_0, Psi_0_inv, nu_0) {
  kappa_q <- kappa_0 + N
  nu_q <- nu_0 + N
  m_q <- (kappa_0 * m_0 + M_sum) / kappa_q

  Psi_q_inv <- Psi_0_inv + B_sum + kappa_0 * tcrossprod(m_0) - kappa_q * tcrossprod(m_q)
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