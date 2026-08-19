.afmvn_components_alpha <- function(data_ctx, alpha, eps_alpha = 1e-8) {
  X <- data_ctx$X
  logX <- data_ctx$logX
  logX_sum <- data_ctx$logX_sum
  D <- data_ctx$D
  N <- data_ctx$N
  Ht <- t(data_ctx$H)
  
  if (abs(alpha) < eps_alpha) {
    W0 <- sweep(logX, 1, rowMeans(logX), "-")
    z0 <- W0 %*% Ht
    logJ0 <- -logX_sum - 0.5 * log(D)
    return(list(
      alpha = alpha,
      near_zero = TRUE,
      z0 = z0,
      z1 = NULL,
      logJ0 = logJ0,
      logJ1 = rep(-Inf, N),
      W_alpha_star = rep(NA, N)
    ))
  }


  A <- alpha*logX
  amax <- apply(A, 1, max)
  log_mean <- amax + log(rowMeans(exp(A - amax)))

  U <- sweep(A, 1, log_mean, "-")
  W0 <- (exp(U) - 1)/alpha

  W_alpha_star <- apply(alpha * W0, 1, min)

  z0 <- W0 %*% Ht
  z1 <- (W0/W_alpha_star^2) %*% Ht

  log_sum <- log_mean + log(D)

  logJ0 <- (D - 0.5) * log(D) + (alpha - 1) * logX_sum - D * log_sum
  logJ1 <- logJ0 - 2 * (D - 1) * log(abs(W_alpha_star))

  list(
    alpha = alpha,
    near_zero = FALSE,
    z0 = z0,
    z1 = z1,
    logJ0 = logJ0,
    logJ1 = logJ1,
    W_alpha_star = W_alpha_star
  )
}

.afmvn_components_eta <- function(data_ctx, eta, eps_alpha = 1e-8) {
  alpha <- tanh(eta)
  comp <- .afmvn_components_alpha(
    data_ctx = data_ctx,
    alpha = alpha,
    eps_alpha = eps_alpha
  )

  comp$eta <- eta
  comp
}


.afmvn_r_stats <- function(data_ctx, gh_ctx, m_eta, s2_eta,
                           quad_ctx, eps_alpha = 1e-8) {
  N <- data_ctx$N

  q <- .eta_quadrature(m_eta = m_eta, s2_eta = s2_eta, gh_ctx = gh_ctx)

  J <- matrix(0, nrow = N, ncol = 2)
  quad <- matrix(0, nrow = N, ncol = 2)

  for (ell in seq_along(q$eta)) {
    comp <- .afmvn_components_eta(data_ctx = data_ctx, eta = q$eta[ell], eps_alpha = eps_alpha)
    w <- q$weight[ell]

    J[, 1] <- J[, 1] + w * comp$logJ0
    quad[, 1] <- quad[, 1] + w * .afmvn_quad_rows(comp$z0, quad_ctx)

    if (!isTRUE(comp$near_zero)) {
      J[, 2] <- J[, 2] + w * comp$logJ1
      quad[, 2] <- quad[, 2] + w * .afmvn_quad_rows(comp$z1, quad_ctx)
    }
  }

  list(J = J, quad = quad)
}