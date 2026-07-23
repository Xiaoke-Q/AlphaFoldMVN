.E_logdet_lambda <- function(Psi_q, nu_q) {
  Psi_q <- as.matrix(Psi_q)
  p <- nrow(Psi_q)
  sum(digamma((nu_q + 1 - seq_len(p)) / 2)) + p * log(2) + as.numeric(determinant(Psi_q, logarithm = TRUE)$modulus)
}

.afmvn_eta_prior_KL <- function(m_eta, s2_eta, a0, s0_sq) {
  -0.5*(log(s0_sq / s2_eta) + (s2_eta + (m_eta - a0)^2) / s0_sq - 1)
}

.afmvn_eta_objective <- function(X, m_eta, s2_eta, r,
                                 m_q, Psi_q, nu_q,
                                 a0 = 0, s0_sq = 10,
                                 n_quad = 20L, eps_alpha = 1e-8) {
  X <- as.matrix(X)
  r <- as.matrix(r)
  Psi_q <- as.matrix(Psi_q)

  N <- nrow(X)
  D <- ncol(X)
  p <- D - 1

  mbj <- .afmvn_components_qeta(X = X, m_eta = m_eta, s2_eta = s2_eta, n_quad = n_quad, eps_alpha = eps_alpha)
  M <- mbj$M
  B <- mbj$B
  J <- mbj$J

  KL_term <- .afmvn_eta_prior_KL(m_eta = m_eta, s2_eta = s2_eta, a0 = a0, s0_sq = s0_sq)
  j_term <- sum(r*J)

  quad_term <- 0
  for (i in 1:N) {
    for (g in 1:2) {
      M_ig <- M[i, g, ]
      B_ig <- B[i, g, , ]

      quad_ig <- sum(Psi_q * t(B_ig)) -
        2*as.numeric(t(m_q)%*%Psi_q%*%M_ig) +
        as.numeric(t(m_q)%*%Psi_q%*%m_q)
      quad_term <- quad_term + r[i, g] * quad_ig
    }
  }

  KL_term + j_term - 0.5*nu_q*quad_term
}



.afmvn_update_eta <- function(X, m_eta, s2_eta, r,
                              m_q, Psi_q, nu_q,
                              a0 = 0, s0_sq = 10,
                              n_quad = 20L,
                              eps_alpha = 1e-8,
                              control = list(maxit = 300L)) {
  objective_neg <- function(par){
    m_eta_cur <- par[1L]
    s2_eta_cur <- exp(par[2L])

    value <- .afmvn_eta_objective(
      X = X, m_eta = m_eta_cur, s2_eta = s2_eta_cur, r = r,
      m_q = m_q, Psi_q = Psi_q, nu_q = nu_q,
      a0 = a0, s0_sq = s0_sq,
      n_quad = n_quad,
      eps_alpha = eps_alpha
    )
    -value
  }

  par0 <- c(m_eta, log(s2_eta))
  value_old <- -objective_neg(par0)

  opt <- stats::optim(
    par = par0,
    fn = objective_neg,
    method = "Nelder-Mead",
    control = control
  )

  m_eta_new <- opt$par[1]
  s2_eta_new <- exp(opt$par[2])
  value_new <- -opt$value

  list(
    m_eta = m_eta_new,
    s2_eta = s2_eta_new,
    log_s2_eta = opt$par[2],
    objective_old = value_old,
    objective = value_new
  )
}

