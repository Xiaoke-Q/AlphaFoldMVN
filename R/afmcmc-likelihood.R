loglik_af <- function(X, alpha, mu, Sigma) {
  ll <- dafmvn(X, alpha = alpha, mean = mu, sigma = Sigma, log = TRUE)
  if (length(ll) > 1) sum(ll) else as.numeric(sum(ll))
}