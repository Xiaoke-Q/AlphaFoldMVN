alpha_to_z <- function(alpha, eps = 1e-6){
  u <- (alpha + 1) / 2
  u <- pmin(pmax(u, eps), 1 - eps)
  log(u) - log(1 - u)
}

z_to_alpha <- function(z) {
  u <- 1 / (1 + exp(-z))
  2 * u - 1
}

logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

reflect_to_interval <- function(x) {
  width <- 2
  y <- (x + 1)%%(2*width)
  if (y > width){
    y <- 2*width - y
  }
    
  y-1 
}

Sigma_to_params <- function(Sigma){
  d <- nrow(Sigma)
  L <- t(chol(Sigma))
  params <- numeric(0)
  for (i in 1:d) {
    if (i > 1) params <- c(params, L[i, 1:(i-1)])
    params <- c(params, log(L[i, i]))
  }

  params
}

params_to_Sigma <- function(params, d){
  L <- matrix(0, d, d)
  k <- 1
  for (i in 1:d){
    if (i > 1){
      L[i, 1:(i-1)] <- params[k:(k + i - 2)]
      k <- k + i - 1
    }
    L[i, i] <- exp(params[k])
    k <- k + 1
  }
  L%*%t(L)
}