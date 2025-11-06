#' The alpha-folding multivariate normal distribution
#'
#' @param n A scalar. Number of observations.
#' @param alpha A numeric parameter of the \eqn{\alpha}-transformation applied to compositional data.
#' @param mean A vector. The mean of latent multivariate normal distribution.
#' @param sigma A symmetric psitive-definite matrix.
#'
#' @returns A matrix of n rows. Each row is an observation. 
#' @name afmvn
#' @export
#' @examples rafmvn(n = 50, alpha = 0.7, mean = c(-1,1), sigma = diag(1, 2))
rafmvn <- function(n = 1, alpha, mean, sigma){
  Y <- mvtnorm::rmvnorm(n = n, mean = mean, sigma = sigma)
  X <- alpha_fold_inverse(Y, alpha = alpha)
  X
}

#' Density function of alpha-folding multivariate normal distribution.
#'
#' @param X A n-by-D matrix.
#' @param alpha A numeric parameter of alpha-folding multivariate normal distribution.
#' @param mean A vector. The mean of latent multivariate normal distribution.
#' @param sigma A symmetric psitive-definite matrix.
#' @param log Logical; return log-likelihood function when log = TRUE.
#'
#' @returns A vector of density values for every observation. 
#'
#' @export
#' @examples dafmvn(X = X, alpha = 0.7, mean = rep(1, D-1), sigma = diag(1, D-1))
dafmvn <- function(X, alpha, mean, sigma, log = TRUE){
  D <- ncol(X)
  N <- nrow(X)
  W <- alpha_trans(X, alpha)
  W_alpha_star <- apply(W*alpha, MARGIN = 1, min) 
  z_alpha_1 <- alpha_fold(X = X, alpha = alpha, unfold = FALSE)
  z_alpha_2 <- alpha_fold(X = X, alpha = alpha, unfold = TRUE)

  if (length(mean) != (D-1))
    stop("length(mean) must equal the dimension of X.")
  if (!is.matrix(sigma) || any(dim(sigma) != c(D-1, D-1)))
    stop("sigma must be a (D-1)-by-(D-1) symmetric positive definite matrix.")
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix.")
  }
  const_term <- -0.5*(D-1)*log(2*pi)
  jacobian_term <- (D-1+0.5)*log(D) + (alpha-1)*rowSums(log(X)) - D*log(rowSums(X^alpha))
  sig_term <- -0.5 * as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
  A <- abs(W_alpha_star)^(-2*(D-1))
  Q_1 <- 0.5*mahalanobis(z_alpha_1, center = mean, cov = sigma)
  Q_2 <- 0.5*mahalanobis(z_alpha_2, center = mean, cov = sigma)
  if(alpha == 0){
    loglik <- const_term +sig_term - Q_1 - rowSums(log(X)) - 0.5*log(D)
  }else{
    loglik <- const_term + jacobian_term + sig_term - Q_1 + log(1+A*exp(Q_1 - Q_2))
  }
  if(isTRUE(log)){
    loglik
  }else(
    exp(loglik)
  )
}