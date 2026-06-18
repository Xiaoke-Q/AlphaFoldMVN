# Construct Gauss-Hermite quadrature nodes and weights for q(eta) = N(m_eta, s2_eta).
.eta_quadrature <- function(m_eta, s2_eta, n_quad = 20L) {
  if (!is.numeric(s2_eta) || length(s2_eta) != 1L || !is.finite(s2_eta) || s2_eta <= 0) {
    stop("s2_eta must be a positive finite numeric value", call. = FALSE)
  }

  gh <- fastGHQuad::gaussHermiteData(as.integer(n_quad))

  eta <- m_eta + sqrt(2 * s2_eta) * gh$x
  weight <- gh$w / sqrt(pi)

  list(
    eta = eta,
    weight = weight,
    raw_nodes = gh$x,
    raw_weights = gh$w
  )
}

# Approximate E[f(eta)] under q(eta) = N(m_eta, s2_eta).
.gh_expectation <- function(f, m_eta, s2_eta, n_quad = 20L) {
  q <- .eta_quadrature(m_eta, s2_eta, n_quad = n_quad)
  vals <- vapply(q$eta, f, numeric(1L))
  sum(q$weight * vals)
}