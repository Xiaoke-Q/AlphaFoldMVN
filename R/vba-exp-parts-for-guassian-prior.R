.afmvn_components_alpha <- function(data_ctx, alpha, eps_alpha = 1e-8) {
  X <- data_ctx$X
  logX <- data_ctx$logX
  logX_sum <- data_ctx$logX_sum
  D <- data_ctx$D
  N <- data_ctx$N
  Ht <- data_ctx$Ht
  
  if (abs(alpha) < eps_alpha) {
    z0 <- alpha_fold(X = X, alpha = 0, unfold = FALSE)
    logJ0 <- -rowSums(log(X)) - 0.5 * log(D)
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

.afmvn_branch_sum <- function(M, B, J, z, logJ, weight, branch) {

  N <- nrow(z)

  M[, branch, ] <- M[, branch, ] + weight * z
  J[, branch] <- J[, branch] + weight * logJ

  p <- ncol(z)
  for (a in seq_len(p)) {
    for (b in seq_len(p)) {
      B[, branch, a, b] <- B[, branch, a, b] + weight * z[, a] * z[, b]
    }
  }

  list(M = M, B = B, J = J)
}

.afmvn_components_qeta <- function(data_ctx, gh_ctx, m_eta, s2_eta, n_quad = 20L, eps_alpha = 1e-8) {
  N <- data_ctx$N
  p <- data_ctx$p

  q <- .eta_quadrature(
    m_eta = m_eta,
    s2_eta = s2_eta,
    gh_ctx = gh_ctx
  )

  M <- array(0, dim = c(N, 2, p))
  B <- array(0, dim = c(N, 2, p, p))
  J <- matrix(0, nrow = N, ncol = 2)

  near_zero_node <- logical(length(q$eta))

  for (ell in seq_along(q$eta)) {
    comp <- .afmvn_components_eta(data_ctx = data_ctx, eta = q$eta[ell], eps_alpha = eps_alpha)

    near_zero_node[ell] <- isTRUE(comp$near_zero)
    w <- q$weight[ell]

    acc_sum <- .afmvn_branch_sum(M = M, B = B,J = J, z = comp$z0, logJ = comp$logJ0, weight = w, branch = 1)

    M <- acc_sum$M
    B <- acc_sum$B
    J <- acc_sum$J

    if (!isTRUE(comp$near_zero)) {
      acc_sum_2 <- .afmvn_branch_sum(M = M, B = B,J = J, z = comp$z1, logJ = comp$logJ1, weight = w, branch = 2)
      M <- acc_sum_2$M
      B <- acc_sum_2$B
      J <- acc_sum_2$J
    }
  }

  list(
    M = M,
    B = B,
    J = J,
    eta = q$eta,
    alpha = tanh(q$eta),
    weight = q$weight,
    near_zero_node = near_zero_node
  )
}