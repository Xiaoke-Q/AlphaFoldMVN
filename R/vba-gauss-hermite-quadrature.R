# Construct Gauss-Hermite quadrature nodes and weights for q(eta) = N(m_eta, s2_eta).
.eta_quadrature <- function(m_eta, s2_eta, gh_ctx) {
  if (!is.numeric(s2_eta) || length(s2_eta) != 1L || !is.finite(s2_eta) || s2_eta <= 0) {
    stop("s2_eta must be a positive finite numeric value", call. = FALSE)
  }

  list(
    eta = m_eta + sqrt(2 * s2_eta) * gh_ctx$x,
    weight = gh_ctx$weight
  )
}