.afmvn_vb_step <- function(data_ctx, gh_ctx,
                           m_eta, s2_eta,
                           m_q, Psi_q, nu_q,
                           m_0, kappa_0, Psi_0_inv, nu_0,
                           a_0 = 0, s0_sq = 10,
                           eps_alpha = 1e-8,
                           eta_control = list(maxit = 300L)) {
  r_quad_ctx <- .afmvn_quad_context(m_q = m_q, Psi_q = Psi_q)

  r_stats <- .afmvn_r_stats(
    data_ctx = data_ctx,
    gh_ctx = gh_ctx,
    m_eta = m_eta,
    s2_eta = s2_eta,
    quad_ctx = r_quad_ctx,
    eps_alpha = eps_alpha
  )  

  r_out <- .afmvn_update_r(
    J = r_stats$J,
    quad = r_stats$quad,
    nu_q = nu_q
  )  

  nw_stats <- .afmvn_nw_stats(
    data_ctx = data_ctx,
    gh_ctx = gh_ctx,
    m_eta = m_eta,
    s2_eta = s2_eta,
    r = r_out$r,
    eps_alpha = eps_alpha
  )  

  nw_out <- .afmvn_update_nw(
    M_sum = nw_stats$M_sum,
    B_sum = nw_stats$B_sum,
    N = data_ctx$N,
    m_0 = m_0,
    kappa_0 = kappa_0,
    Psi_0_inv = Psi_0_inv,
    nu_0 = nu_0
  )
  
  eta_quad_ctx <- .afmvn_quad_context(
    m_q = nw_out$m_q,
    Psi_q = nw_out$Psi_q
  )

  eta_out <- .afmvn_update_eta(
    data_ctx = data_ctx,
    gh_ctx = gh_ctx,
    m_eta = m_eta,
    s2_eta = s2_eta,
    r = r_out$r,
    quad_ctx = eta_quad_ctx,
    nu_q = nw_out$nu_q,
    a0 = a_0,
    s0_sq = s0_sq,
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

#' Variational Bayes for the AFMVN model
#'
#' Fits the AFMVN model using variational Bayes.
#'
#' @param X Data matrix.
#' @param m_eta_init Initial mean of q(eta).
#' @param s2_eta_init Initial variance of q(eta).
#' @param m_q_init Initial mean vector of q(mu, Lambda).
#' @param Psi_q_init Initial Wishart scale matrix.
#' @param nu_q_init Initial Wishart degrees of freedom.
#' @param m_0 Prior mean vector.
#' @param kappa_0 Prior scaling parameter.
#' @param Psi_0 Prior Wishart scale matrix.
#' @param nu_0 Prior Wishart degrees of freedom.
#' @param a_0 Prior mean of eta.
#' @param s0_sq Prior variance of eta.
#' @param n_quad Number of quadrature nodes.
#' @param eps_alpha Small value used for numerical stability.
#' @param eta_control Control list passed to the eta optimization.
#' @param max_iter Maximum number of VB iterations.
#' @param tol Convergence tolerance.
#'
#' @return A list containing the variational parameter estimates, iteration trace,
#' number of iterations, and convergence indicator.
#'
#' @export
afmvn_vba <- function(X,
                      m_eta_init, s2_eta_init,
                      m_q_init, Psi_q_init, nu_q_init,
                      m_0, kappa_0, Psi_0, nu_0,
                      a_0 = 0, s0_sq = 10,
                      n_quad = 20L,
                      eps_alpha = 1e-8,
                      eta_control = list(maxit = 300L),
                      max_iter = 150, tol = 1e-3) {
  X <- as.matrix(X)
  if (any(X <= 0)) stop("X must be strictly positive.")
  if (any(abs(rowSums(X) - 1) > 1e-8)) stop("Rows of X must sum to 1.")
  D <- ncol(X)
  p <- D - 1
  
  logX <- log(X)
  data_ctx <- list(
    X = as.matrix(X),
    logX = logX,
    logX_sum = rowSums(logX),
    H = helmert(D),
    D = D,
    N = nrow(X),
    p = p
  )
  
  gh <- fastGHQuad::gaussHermiteData(n_quad)
  gh_ctx <- list(
    x = gh$x,
    weight = gh$w/sqrt(pi)
  )
  
  Psi_0_inv <- solve(Psi_0)


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
      data_ctx = data_ctx,
      gh_ctx = gh_ctx,
      m_eta = state_old$m_eta,
      s2_eta = state_old$s2_eta,
      m_q = state_old$m_q,
      Psi_q = state_old$Psi_q,
      nu_q = state_old$nu_q,
      m_0 = m_0,
      kappa_0 = kappa_0,
      Psi_0_inv = Psi_0_inv,
      nu_0 = nu_0,
      a_0 = a_0,
      s0_sq = s0_sq,
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