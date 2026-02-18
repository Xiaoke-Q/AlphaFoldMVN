#' MLE of alpha-folding multivariate normal distribution given alpha
#'
#' @param X A n-by-D matrix.
#' @param alpha A numeric parameter of alpha-folding multivariate normal distribution.
#' @param tol A threshold of alpha. When |alpha| < tol, the clt is adopted.
#' @param max.iter The number of alternating iterations of solving for mu and sigma. 
#'
#' @returns A list of estimated mu and sigma.
#'
#' @export
#' @examples afmvn_mle(X = X, alpha = 0.7, tol = 1e-6, max.iter = 150)
afmvn_mle <- function(X, 
                      alpha, 
                      tol = 1e-6, 
                      max.iter = 150){
  if(abs(alpha) < tol){
    D <- ncol(X)
    N <- nrow(X)
    H <- helmert(D)
    W <- alpha_trans(X, alpha)
    Y0 <- W %*% t(H)
    
    mu_hat <- colMeans(Y0)
    Sigma_hat <- cov(Y0)

    return(list(Sigma_hat = Sigma_hat, mu_hat = mu_hat))
  }

  k_factor <- function(mu, sigma){
    A <- abs(W_alpha_star)^(-2*(D-1))
    Q_1 <- 0.5*mahalanobis(z_alpha_1, center = mu, cov = sigma)
    Q_2 <- 0.5*mahalanobis(z_alpha_2, center = mu, cov = sigma)
    A*exp(Q_1 - Q_2)
  }

  D <- ncol(X)
  N <- nrow(X)
  H <- helmert(D)
  W <- alpha_trans(X, alpha)
  W_alpha_star <- apply(W * alpha, MARGIN = 1, min) 
  z_alpha_1 <- W %*% t(H)
  z_alpha_2 <- z_alpha_1/W_alpha_star^2
  
  
  mu_hat <- (colMeans(z_alpha_2) + colMeans(z_alpha_2))/2
  Sigma_hat <- (cov(z_alpha_1) + cov(z_alpha_2))/2

  count <- 0
  flag = TRUE
  while(count <= max.iter){
    K <- k_factor(mu = mu_hat, sigma = Sigma_hat)
    
    w1 <- 1/(1+K)
    w2 <- 1- w1
    mu_hat_updated <- colMeans(w1 * z_alpha_1 + w2 * z_alpha_2)
    z_alpha_1_shifted_weighted <- sweep(z_alpha_1, 2, mu_hat_updated, "-")*sqrt(w1)
    z_alpha_2_shifted_weighted <- sweep(z_alpha_2, 2, mu_hat_updated, "-")*sqrt(w2)
    Sigma_hat_updated <- (crossprod(z_alpha_1_shifted_weighted) + crossprod(z_alpha_2_shifted_weighted))/N
    mu_hat <- mu_hat_updated
    Sigma_hat <- Sigma_hat_updated
    count = count + 1
  }
  list(Sigma_hat = Sigma_hat, mu_hat = mu_hat)
}


