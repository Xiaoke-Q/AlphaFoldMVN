.afmvn_vb_step <- function(X,
                           m_eta, s2_eta,
                           m_q, kappa_q, Psi_q, nu_q,
                           m_0, kappa_0, Psi_0, nu_0,
                           a_0 = 0, s0_sq = 10,
                           n_quad = 20L,
                           eps_alpha = 1e-8,
                           eta_control = list(maxit = 300L)) {
  mbj <- .afmvn_components_qeta(
    X = X,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = n_quad,
    eps_alpha = eps_alpha
  )

  r_out <- .afmvn_update_r(
    M = mbj$M,
    B = mbj$B,
    J = mbj$J,
    m_q = m_q,
    Psi_q = Psi_q,
    nu_q = nu_q
  )

  nw_out <- .afmvn_update_nw(
    M = mbj$M,
    B = mbj$B,
    r = r_out$r,
    m_0 = m_0,
    kappa_0 = kappa_0,
    Psi_0 = Psi_0,
    nu_0 = nu_0
  )

  eta_out <- .afmvn_update_eta(
    X = X,
    m_eta = m_eta,
    s2_eta = s2_eta,
    r = r_out$r,
    m_q = nw_out$m_q,
    Psi_q = nw_out$Psi_q,
    nu_q = nw_out$nu_q,
    a0 = a_0,
    s0_sq = s0_sq,
    n_quad = n_quad,
    eps_alpha = eps_alpha,
    control = eta_control
  )

  list(
    m_eta = eta_out$m_eta,
    s2_eta = eta_out$s2_eta,
    r = r_out$r,
    m_q = nw_out$m_q,
    kappa_q = nw_out$kappa_q,
    Psi_q = nw_out$Psi_q,
    nu_q = nw_out$nu_q,
    eta_objective = eta_out$objective,
    eta_convergence = eta_out$convergence
  )
}


.afmvn_relative_change <- function(new, old) {
  if (is.null(old)) {
    return(Inf)
  }

  max(abs(new - old) / (1 + abs(old)))
}

afmvn_vba <- function(X,
                      m_eta_init, s2_eta_init,
                      m_q_init, Psi_q_init, nu_q_init,
                      m_0, kappa_0, Psi_0, nu_0,
                      a_0 = 0, s0_sq = 10,
                      n_quad = 20L,
                      eps_alpha = 1e-8,
                      eta_control = list(maxit = 300L),
                      max_iter = 150, tol = 1e-3) {
  state <- list(
    m_eta = m_eta_init,
    s2_eta = s2_eta_init,
    r = NULL,
    m_q = m_q_init,
    kappa_q = kappa_0,
    Psi_q = Psi_q_init,
    nu_q = nu_q_init
  )

  trace <- vector("list", max_iter)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    state_old <- state
    state <- .afmvn_vb_step(
      X = X,
      m_eta = state_old$m_eta,
      s2_eta = state_old$s2_eta,
      m_q = state_old$m_q,
      Psi_q = state_old$Psi_q,
      nu_q = state_old$nu_q,
      m_0 = m_0,
      kappa_0 = kappa_0,
      Psi_0 = Psi_0,
      nu_0 = nu_0,
      a_0 = a_0,
      s0_sq = s0_sq,
      n_quad = n_quad,
      eps_alpha = eps_alpha,
      eta_control = eta_control
    )


    delta_eta <- .afmvn_relative_change(
      c(state$m_eta, log(state$s2_eta)),
      c(state_old$m_eta, log(state_old$s2_eta))
    )

    delta_nw <- .afmvn_relative_change(
      c(state$m_q, as.vector(state$Psi_q), state$nu_q),
      c(state_old$m_q, as.vector(state_old$Psi_q), state_old$nu_q)
    )

    delta_r <- .afmvn_relative_change(
      state$r,
      state_old$r
    )

    max_change <- max(delta_eta, delta_nw, delta_r)

    trace[[iter]] <- list(
      iter = iter,
      m_eta = state$m_eta,
      alpha = tanh(state$m_eta),
      s2_eta = state$s2_eta,
      eta_objective = state$eta_objective,
      delta_eta = delta_eta,
      delta_nw = delta_nw,
      delta_r = delta_r,
      max_change = max_change
    )
    cat(sprintf("iter = %d, alpha = %.6f,      s2_eta = %.6e, max_change = %.6e\n",
                iter,       tanh(state$m_eta), state$s2_eta,  max_change))
    if (iter > 1L && max_change < tol) {
      converged <- TRUE
      break
    }
  }

  trace <- do.call(rbind, trace[seq_len(iter)])

  list(
    m_eta = state$m_eta,
    s2_eta = state$s2_eta,
    r = state$r,
    m_q = state$m_q,
    kappa_q = state$kappa_q,
    Psi_q = state$Psi_q,
    nu_q = state$nu_q,
    trace = trace,
    iterations = iter,
    converged = converged
  )
}