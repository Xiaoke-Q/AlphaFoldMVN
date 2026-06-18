
test_that(".eta_quadrature validates s2_eta", {
  expect_error(
    .eta_quadrature(m_eta = 0, s2_eta = 0),
    "s2_eta must be a positive finite numeric scalar"
  )

  expect_error(
    .eta_quadrature(m_eta = 0, s2_eta = -1),
    "s2_eta must be a positive finite numeric scalar"
  )

  expect_error(
    .eta_quadrature(m_eta = 0, s2_eta = Inf),
    "s2_eta must be a positive finite numeric scalar"
  )

  expect_error(
    .eta_quadrature(m_eta = 0, s2_eta = c(1, 2)),
    "s2_eta must be a positive finite numeric scalar"
  )
})

test_that(".gh_expectation integrates constants", {
  val <- .gh_expectation(
    f = function(eta) 3,
    m_eta = 0.5,
    s2_eta = 2,
    n_quad = 20L
  )

  expect_equal(val, 3, tolerance = 1e-12)
})

test_that(".gh_expectation matches normal first and second moments", {
  m_eta <- 1.3
  s2_eta <- 0.8

  e1 <- .gh_expectation(
    f = function(eta) eta,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 20L
  )

  e2 <- .gh_expectation(
    f = function(eta) eta^2,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 20L
  )

  expect_equal(e1, m_eta, tolerance = 1e-12)
  expect_equal(e2, s2_eta + m_eta^2, tolerance = 1e-12)
})

test_that(".gh_expectation matches normal fourth moment", {
  m_eta <- -0.4
  s2_eta <- 1.7

  e4 <- .gh_expectation(
    f = function(eta) eta^4,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 20L
  )

  target <- m_eta^4 + 6 * m_eta^2 * s2_eta + 3 * s2_eta^2

  expect_equal(e4, target, tolerance = 1e-10)
})

test_that(".gh_expectation handles nonlinear smooth functions", {
  m_eta <- 0.2
  s2_eta <- 0.5

  val <- .gh_expectation(
    f = function(eta) exp(eta),
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 30L
  )

  target <- exp(m_eta + 0.5 * s2_eta)

  expect_equal(val, target, tolerance = 1e-10)
})

test_that(".gh_expectation for E[tanh(eta)] is stable as n_quad increases", {
  m_eta <- 0.7
  s2_eta <- 0.4

  val20 <- .gh_expectation(
    f = tanh,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 20L
  )

  val40 <- .gh_expectation(
    f = tanh,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 40L
  )

  val80 <- .gh_expectation(
    f = tanh,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = 80L
  )

  expect_equal(val20, val40, tolerance = 1e-8)
  expect_equal(val40, val80, tolerance = 1e-8)
})

test_that(".gh_expectation respects odd symmetry of tanh", {
  val <- .gh_expectation(
    f = tanh,
    m_eta = 0,
    s2_eta = 2,
    n_quad = 30L
  )

  expect_equal(val, 0, tolerance = 1e-14)
})