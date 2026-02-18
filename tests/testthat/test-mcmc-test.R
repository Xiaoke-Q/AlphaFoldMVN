testthat::test_that("MCMC test runs and returns correct shapes", {
  testthat::skip_if_not_installed("mvtnorm")
  set.seed(1)

  d <- 2
  n <- 8

  prior <- list(
    m0 = rep(1, d),
    kappa0 = 1,
    Psi = diag(d),
    nu = 6
  )

  init <- list(
    alpha = 0.2,
    mu = rep(1, d),
    Sigma = diag(d)
  )

  # Generate data from the model itself -> guaranteed to be in-domain
  X <- rafmvn(n = n, alpha = init$alpha, mean = init$mu, sigma = init$Sigma)$X

  testthat::expect_equal(dim(X), c(n, d+1))
  testthat::expect_true(all(is.finite(X)))

  # If your model requires strict positivity, keep this:
  testthat::expect_true(all(X > 0))

  # Likelihood finite at init
  ll0 <- dafmvn(X, alpha = init$alpha, mean = init$mu, sigma = init$Sigma, log = TRUE)
  testthat::expect_true(is.finite(sum(ll0)))

  out <- afmvn_mcmc(
    X = X, prior = prior, init = init,
    warm_up = 20, warm_band = 10, n_iter = 20,
    alpha_repropose = 2,
    alpha_step = 0.2, mu_step = 0.1, Sigma_step = 0.05
  )

  testthat::expect_true(is.list(out))
  testthat::expect_true(all(c("alpha","mu","Sigma","tuned","acceptance") %in% names(out)))

  testthat::expect_length(out$alpha, 20)
  testthat::expect_equal(dim(out$mu), c(20, d))
  testthat::expect_equal(dim(out$Sigma), c(d, d, 20))

  for (t in 1:20) {
    testthat::expect_silent(chol(out$Sigma[, , t]))
  }

  testthat::expect_true(all(out$acceptance >= 0))
  testthat::expect_true(all(out$acceptance <= 1))
})