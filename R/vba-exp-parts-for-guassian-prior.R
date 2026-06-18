.afmvn_components_alpha <- function(X, alpha, eps_alpha = 1e-6) {
  D <- ncol(X)
  N <- nrow(X)

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

  W0 <- alpha_trans(X, alpha = alpha, unfold = FALSE)
  W_alpha_star <- apply(alpha*W0, MARGIN = 1, min)

  z0 <- alpha_fold(X = X, alpha = alpha, unfold = FALSE)
  z1 <- alpha_fold(X = X, alpha = alpha, unfold = TRUE)


  logJ0 <- (D - 0.5)*log(D) + (alpha - 1)*rowSums(log(X)) - D*log(rowSums(X^alpha))
  logJ1 <- logJ0 - 2*(D - 1) * log(abs(W_alpha_star))

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

.afmvn_components_eta <- function(X, eta, eps_alpha = 1e-8) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta)) {
    stop("eta must be a finite numeric scalar.", call. = FALSE)
  }

  alpha <- tanh(eta)
  comp <- .afmvn_components_alpha(
    X = X,
    alpha = alpha,
    eps_alpha = eps_alpha
  )

  comp$eta <- eta
  comp
}

.afmvn_branch_sum <- function(M, B, J, z, logJ, weight, branch) {

  N <- nrow(z)

  J[, branch] <- J[, branch] + weight*logJ

  for (i in 1:N) {
    z_i <- z[i,]
    M[i, branch, ] <- M[i, branch, ] + weight*z_i
    B[i, branch, , ] <- B[i, branch, , ] + weight*tcrossprod(z_i)
  }

  list(M = M, B = B, J = J)
}

.afmvn_components_qeta <- function(X, m_eta, s2_eta, n_quad = 20L, eps_alpha = 1e-8) {
  X <- as.matrix(X)
  N <- nrow(X)
  D <- ncol(X)
  p <- D - 1

  q <- .eta_quadrature(
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = n_quad
  )

  M <- array(0, dim = c(N, 2, p))
  B <- array(0, dim = c(N, 2, p, p))
  J <- matrix(0, nrow = N, ncol = 2)

  near_zero_node <- logical(length(q$eta))

  for (ell in seq_along(q$eta)) {
    comp <- .afmvn_components_eta(X = X, eta = q$eta[ell], eps_alpha = eps_alpha)

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