propose_alpha <- function(m, s) {
  z_prop <- rnorm(1, mean = m, sd = s)
  z_to_alpha(z_prop)
}

propose_mu <- function(mu, step_sd, Sigma){
  as.numeric(MASS::mvrnorm(1, mu, step_sd^2 * Sigma))
}

propose_Sigma <- function(Sigma, step_sd){
  d <- nrow(Sigma)
  p <- Sigma_to_params(Sigma)
  p_prop <- p + rnorm(length(p), 0, step_sd)
  params_to_Sigma(p_prop, d)
}