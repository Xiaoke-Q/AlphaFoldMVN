test_that(".afmvn_components_alpha returns valid components away from zero", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  alpha <- 0.7
  comp <- .afmvn_components_alpha(X = X, alpha = alpha, eps_alpha = 1e-8)

  expect_false(comp$near_zero)
  expect_equal(comp$alpha, alpha)

  expect_true(is.matrix(comp$z0))
  expect_true(is.matrix(comp$z1))
  expect_equal(dim(comp$z0), c(nrow(X), ncol(X) - 1L))
  expect_equal(dim(comp$z1), c(nrow(X), ncol(X) - 1L))

  expect_length(comp$logJ0, nrow(X))
  expect_length(comp$logJ1, nrow(X))
  expect_length(comp$W_alpha_star, nrow(X))

  expect_true(all(is.finite(comp$z0)))
  expect_true(all(is.finite(comp$z1)))
  expect_true(all(is.finite(comp$logJ0)))
  expect_true(all(is.finite(comp$logJ1)))
  expect_true(all(is.finite(comp$W_alpha_star)))
})

test_that(".afmvn_components_alpha handles near-zero alpha by ordinary branch only", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  comp <- .afmvn_components_alpha(X = X, alpha = 1e-10, eps_alpha = 1e-8)

  expect_true(comp$near_zero)
  expect_true(is.matrix(comp$z0))
  expect_null(comp$z1)

  expect_equal(dim(comp$z0), c(nrow(X), ncol(X) - 1L))
  expect_length(comp$logJ0, nrow(X))
  expect_length(comp$logJ1, nrow(X))
  expect_length(comp$W_alpha_star, nrow(X))

  expect_true(all(is.finite(comp$z0)))
  expect_true(all(is.finite(comp$logJ0)))
  expect_true(all(is.infinite(comp$logJ1)))
  expect_true(all(comp$logJ1 < 0))
  expect_true(all(is.na(comp$W_alpha_star)))
})



test_that(".afmvn_components_alpha computes log Jacobians consistently", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  D <- ncol(X)

  alpha_values <- c(-0.8, -0.4, -0.05, -0.005, 0.005,0.05, 0.4, 0.8)

  for (alpha in alpha_values) {
    comp <- .afmvn_components_alpha(
      X = X,
      alpha = alpha,
      eps_alpha = 1e-8
    )

    W0 <- alpha_trans(X, alpha = alpha, unfold = FALSE)
    W_alpha_star <- apply(alpha * W0, MARGIN = 1L, min)

    logJ0_target <- (D - 0.5) * log(D) +
      (alpha - 1) * rowSums(log(X)) -
      D * log(rowSums(X^alpha))

    logJ1_target <- logJ0_target -
      2 * (D - 1) * log(abs(W_alpha_star))

    expect_false(comp$near_zero, info = paste("alpha =", alpha))
    expect_equal(comp$logJ0, logJ0_target, tolerance = 1e-12, info = paste("alpha =", alpha))
    expect_equal(comp$logJ1, logJ1_target, tolerance = 1e-12, info = paste("alpha =", alpha))
    expect_equal(comp$W_alpha_star, W_alpha_star, tolerance = 1e-12, info = paste("alpha =", alpha))
  }
})

test_that(".afmvn_components_alpha computes alpha-zero log Jacobian", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  D <- ncol(X)
  comp <- .afmvn_components_alpha(X = X, alpha = 0, eps_alpha = 1e-8)

  logJ0_target <- -rowSums(log(X)) - 0.5*log(D)

  expect_equal(comp$logJ0, logJ0_target, tolerance = 1e-12)
})


test_that(".afmvn_components_eta validates eta", {
  X <- matrix(c(0.2, 0.3, 0.5), nrow = 1L)

  expect_error(.afmvn_components_eta(X = X, eta = Inf), "eta must be a finite numeric scalar")
  expect_error(.afmvn_components_eta(X = X, eta = c(0, 1)), "eta must be a finite numeric scalar")
})

test_that(".afmvn_branch_sum accumulates one branch correctly", {
  N <- 2L
  p <- 2L

  M <- array(0, dim = c(N, 2L, p))
  B <- array(0, dim = c(N, 2L, p, p))
  J <- matrix(0, nrow = N, ncol = 2L)

  z <- rbind(
    c(1, 2),
    c(3, 4)
  )

  logJ <- c(0.5, -0.2)
  weight <- 0.25

  out <- .afmvn_branch_sum(
    M = M,
    B = B,
    J = J,
    z = z,
    logJ = logJ,
    weight = weight,
    branch = 2L
  )

  expect_equal(out$J[, 2L], weight*logJ)
  expect_equal(out$J[, 1L], c(0, 0))

  expect_equal(out$M[1L, 2L, ], weight*z[1L, ])
  expect_equal(out$M[2L, 2L, ], weight*z[2L, ])

  expect_equal(out$B[1L, 2L, , ], weight*tcrossprod(z[1L, ]))
  expect_equal(out$B[2L, 2L, , ], weight*tcrossprod(z[2L, ]))

  expect_equal(out$M[, 1L, ], matrix(0, nrow = N, ncol = p))
})

test_that(".afmvn_components_qeta returns valid M, B, J arrays", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  out <- .afmvn_components_qeta(
    X = X,
    m_eta = 0.5,
    s2_eta = 0.1,
    n_quad = 20L,
    eps_alpha = 1e-8
  )

  N <- nrow(X)
  p <- ncol(X) - 1L

  expect_equal(dim(out$M), c(N, 2L, p))
  expect_equal(dim(out$B), c(N, 2L, p, p))
  expect_equal(dim(out$J), c(N, 2L))

  expect_length(out$eta, 20L)
  expect_length(out$alpha, 20L)
  expect_length(out$weight, 20L)
  expect_length(out$near_zero_node, 20L)

  expect_true(all(is.finite(out$M)))
  expect_true(all(is.finite(out$B)))
  expect_true(all(is.finite(out$J)))
  expect_equal(sum(out$weight), 1, tolerance = 1e-12)
})

test_that(".afmvn_components_qeta approximates fixed eta when s2_eta is very small", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  m_eta <- 0.6
  s2_eta <- 1e-10

  out <- .afmvn_components_qeta(
    X = X,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 20L,
    eps_alpha = 1e-8
  )

  comp <- .afmvn_components_eta(
    X = X,
    eta = m_eta,
    eps_alpha = 1e-8
  )

  for (i in seq_len(nrow(X))) {
    expect_equal(out$M[i, 1L, ], comp$z0[i, ], tolerance = 1e-6)
    expect_equal(out$B[i, 1L, , ], tcrossprod(comp$z0[i, ]), tolerance = 1e-6)
  }

  expect_equal(out$J[, 1L], comp$logJ0, tolerance = 1e-6)

  if (!isTRUE(comp$near_zero)) {
    for (i in seq_len(nrow(X))) {
      expect_equal(out$M[i, 2L, ], comp$z1[i, ], tolerance = 1e-6)
      expect_equal(out$B[i, 2L, , ], tcrossprod(comp$z1[i, ]), tolerance = 1e-6)
    }

    expect_equal(out$J[, 2L], comp$logJ1, tolerance = 1e-6)
  }
})

test_that(".afmvn_components_qeta handles near-zero eta region", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2)
  )

  out <- .afmvn_components_qeta(
    X = X,
    m_eta = 0,
    s2_eta = 1e-16,
    n_quad = 20L,
    eps_alpha = 1e-6
  )

  expect_true(all(out$near_zero_node))

  expect_true(all(is.finite(out$M[, 1L, ])))
  expect_true(all(is.finite(out$B[, 1L, , ])))
  expect_true(all(is.finite(out$J[, 1L])))

  expect_equal(out$M[, 2L, ], matrix(0, nrow = nrow(X), ncol = ncol(X) - 1L))
  expect_true(all(out$B[, 2L, , ] == 0))
  expect_equal(out$J[, 2L], rep(0, nrow(X)))
})