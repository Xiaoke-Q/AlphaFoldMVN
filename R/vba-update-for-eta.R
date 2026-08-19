.E_logdet_lambda <- function(Psi_q, nu_q) {
  Psi_q <- as.matrix(Psi_q)
  p <- nrow(Psi_q)
  sum(digamma((nu_q + 1 - seq_len(p)) / 2)) + p * log(2) + as.numeric(determinant(Psi_q, logarithm = TRUE)$modulus)
}

.afmvn_eta_prior_KL <- function(m_eta, s2_eta, a0, s0_sq) {
  -0.5*(log(s0_sq / s2_eta) + (s2_eta + (m_eta - a0)^2) / s0_sq - 1)
}


.afmvn_eta_objective <- function(X, m_eta, s2_eta, r,
                                 quad_ctx, nu_q,
                                 a0 = 0, s0_sq = 10,
                                 n_quad = 20L, eps_alpha = 1e-8) {
  X <- as.matrix(X)
  r <- as.matrix(r)

  mbj <- .afmvn_components_qeta(
    X = X,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = n_quad,
    eps_alpha = eps_alpha
  )

  quad <- .afmvn_quad_from_mbj(mbj = mbj, quad_ctx = quad_ctx)

  KL_term <- .afmvn_eta_prior_KL(
    m_eta = m_eta,
    s2_eta = s2_eta,
    a0 = a0,
    s0_sq = s0_sq
  )

  j_term <- sum(r * mbj$J)
  quad_term <- sum(r * quad)

  KL_term + j_term - 0.5 * nu_q * quad_term
}


.afmvn_update_eta <- function(X, m_eta, s2_eta, r,
                              m_q = NULL, Psi_q = NULL, nu_q = NULL,
                              quad_ctx = NULL,
                              a0 = 0, s0_sq = 10,
                              n_quad = 20L,
                              eps_alpha = 1e-8,
                              control = list(maxit = 300L)) {
  if (is.null(nu_q)) {
    stop("nu_q must be provided.", call. = FALSE)
  }

  if (is.null(quad_ctx)) {
    if (is.null(m_q) || is.null(Psi_q)) {
      stop("Provide either quad_ctx or both m_q and Psi_q.", call. = FALSE)
    }

    quad_ctx <- .afmvn_quad_context(m_q = m_q, Psi_q = Psi_q)
  }
  
  objective_neg <- function(par){
    m_eta_cur <- par[1L]
    s2_eta_cur <- exp(par[2L])

    value <- .afmvn_eta_objective(
      X = X,
      m_eta = m_eta_cur,
      s2_eta = s2_eta_cur,
      r = r,
      quad_ctx = quad_ctx,
      nu_q = nu_q,
      a0 = a0,
      s0_sq = s0_sq,
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

