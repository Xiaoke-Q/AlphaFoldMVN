r1_prior_diff <- function(alpha_prop, alpha_now, m, s) {
  z_now  <- alpha_to_z(alpha_now)
  z_prop <- alpha_to_z(alpha_prop)

  dnorm(z_now,  mean = m, sd = s, log = TRUE) -
    dnorm(z_prop, mean = m, sd = s, log = TRUE) +
    log1p(-alpha_prop^2) - log1p(-alpha_now^2)
}

r2_prior_diff <- function(mu_prop, mu_now, Sigma_now, prior) {
  Sigma_inv <- chol2inv(chol(Sigma_now))

  dprop <- mu_prop - prior$m0
  dnow  <- mu_now  - prior$m0

  qprop <- mahalanobis(mu_prop, center = prior$m0, cov = Sigma_inv, inverted = TRUE)
  qnow  <- mahalanobis(mu_now,  center = prior$m0, cov = Sigma_inv, inverted = TRUE)

  0.5 * prior$kappa0 * (qnow - qprop)
}
r3_prior_diff <- function(Sigma_prop, Sigma_now, mu_now, prior){
  d <- length(mu_now)
  diff <- mu_now - prior$m0

  cholP <- tryCatch(chol(Sigma_prop), error = function(e) NULL)
  if (is.null(cholP)) return(-Inf)
  cholN <- tryCatch(chol(Sigma_now), error = function(e) NULL)
  if (is.null(cholN)) return(-Inf)

  logdetP <- 2 * sum(log(diag(cholP)))
  logdetN <- 2 * sum(log(diag(cholN)))

  invP <- chol2inv(cholP)
  invN <- chol2inv(cholN)

  trP <- sum(diag(prior$Psi %*% invP))
  trN <- sum(diag(prior$Psi %*% invN))

  uP <- backsolve(cholP, diff, transpose = TRUE)
  qP <- sum(uP^2)
  uN <- backsolve(cholN, diff, transpose = TRUE)
  qN <- sum(uN^2)

  lpP <- -0.5 * (prior$nu + d + 2) * logdetP - 0.5 * trP - 0.5 * prior$kappa0 * qP
  lpN <- -0.5 * (prior$nu + d + 2) * logdetN - 0.5 * trN - 0.5 * prior$kappa0 * qN

  lpP - lpN
}