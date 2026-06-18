test_that(".afmvn_eta_prior_KL is maximized at the prior", {
  a0 <- 0.3
  s0_sq <- 2.5

  val_at_prior <- .afmvn_eta_prior_KL(
    m_eta = a0,
    s2_eta = s0_sq,
    a0 = a0,
    s0_sq = s0_sq
  )

  val_away <- .afmvn_eta_prior_KL(
    m_eta = a0 + 1,
    s2_eta = s0_sq,
    a0 = a0,
    s0_sq = s0_sq
  )

  val_wrong_var <- .afmvn_eta_prior_KL(
    m_eta = a0,
    s2_eta = 0.5 * s0_sq,
    a0 = a0,
    s0_sq = s0_sq
  )

  expect_equal(val_at_prior, 0, tolerance = 1e-12)
  expect_lt(val_away, val_at_prior)
  expect_lt(val_wrong_var, val_at_prior)
})

test_that(".afmvn_eta_objective returns a finite scalar", {
  X <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.1, 0.7, 0.2),
    c(0.4, 0.4, 0.2)
  )

  p <- ncol(X) - 1L

  r <- matrix(c(
    0.8, 0.2,
    0.7, 0.3,
    0.9, 0.1
  ), nrow = nrow(X), byrow = TRUE)

  val <- .afmvn_eta_objective(
    X = X,
    m_eta = 0.4,
    s2_eta = 0.2,
    r = r,
    m_q = rep(0, p),
    Psi_q = diag(p),
    nu_q = p + 5,
    a0 = 0,
    s0_sq = 10,
    n_quad = 20L,
    eps_alpha = 1e-8
  )

  expect_true(is.numeric(val))
  expect_length(val, 1L)
  expect_true(is.finite(val))
})
