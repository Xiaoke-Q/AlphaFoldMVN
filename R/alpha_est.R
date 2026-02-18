#' Estimation of alpha
#'
#' @param X A n-by-D matrix.
#' @param by increment of the alphas
#'
#' @returns list: alpha_hat is the estimated alpha, loglik is a vector of log-likelihood for each alpha 
#'
#' @export
#' @examples alpha_est(X = X, by = 0.01)
alpha_est <- function(X, by = 0.01, level = NULL){
  alphas <- seq(-1, 1, by = by)
  CI = NULL
  loglik <- future.apply::future_sapply(
    alphas,
    FUN = function(a) {
      r <- afmvn_mle(X = X, alpha = a)
      loglik <- sum(dafmvn(X = X, alpha = a,
               mean = r$mu_hat,
               sigma = r$Sigma_hat))
      loglik
    }
  )
  
  m <- which.max(loglik)
  alpha_hat <- alphas[m]
  if(!is.null(level)){
    cut <- max(loglik) - 0.5 * qchisq(level, df = 1)
    diff <- (loglik - cut) > 0
    cross <- c(1, which(diff[2:(length(diff)-1)] != diff[1:(length(diff)-2)]) + 1, length(diff))
    cross_index <- findInterval(m, cross)
    lower_bound <- ifelse(cross_index == m, alphas[m], alphas[cross[cross_index]])
    upper_bound <- ifelse(cross_index == m, alphas[m], alphas[cross[cross_index+1]])
  }
  list(alpha_hat = alpha_hat, loglik = loglik, alphas = alphas, CI = c(lower_bound, upper_bound))
}

#' Estimation of alpha (Method 2)
#'
#' @param X A n-by-D matrix.
#' @param Sigma known sigma
#' @param mu kown mean
#' @returns list: alpha_hat is the estimated alpha, loglik is a vector of log-likelihood for each alpha 
#'
#' @export
#' @examples alpha_est(X = X, by = 0.01)
alpha_est2 <- function(X, mu, sigma, alpha){

  z_alpha_update <- function(alpha){
    z_alpha <<- alpha_fold(X = X, alpha = alpha, unfold = FALSE)
  }
  y_alpha_update <- function(alpha){
    y_alpha <<- alpha_fold(X = X, alpha = alpha, unfold = TRUE)
  }
  w_alpha_update <- function(alpha){
    w_alpha <<- alpha_trans(X, alpha)
  }
  w_alpha_star_update <- function(alpha){
    w_alpha_star <<- apply(w_alpha*alpha, MARGIN = 1, min) 
  }
  k_factor_update <- function(alpha, mu, sigma){
    A <- abs(w_alpha_star)^(-2*(D-1))
    Q_1 <- 0.5*mahalanobis(z_alpha, center = mu, cov = sigma)
    Q_2 <- 0.5*mahalanobis(y_alpha, center = mu, cov = sigma)
    K_alpha <<- A*exp(Q_1 - Q_2)
  }
  u_alpha_update <- function(alpha){
    A <- alpha * logX              
    m <- apply(A, 1L, max)
    shifted_exp <- rowMeans(exp(A - m))
    denom_log <- m + log(shifted_exp) 

    U <- sweep(A, 1L, denom_log, FUN = "-") 
    u_alpha <<- exp(U)
  }
  partial_u_alpha_update <- function(alpha){
    partial_u_alpha <<- (logX - rowSums(u_alpha*logX)) * u_alpha
  }
  A_alpha_update <- function(alpha){

    sr <- t(solve(sigma, t(y_alpha)))
    s <- rowSums((y_alpha - MU) * sr)

    j_star <- apply(u_alpha, MARGIN = 1, which.min)
    u_alpha_star <- apply(u_alpha, MARGIN = 1, min)
    w_part_up <- 2*D*partial_u_alpha[cbind(seq_len(N), j_star)]
    w_part_down <- D*u_alpha_star - 1
    w_part <- w_part_up/w_part_down

    term1 <- sum((1-D*u_alpha)*logX) 
    term2 <- K_alpha/(1+K_alpha) * (D-1 - s) * w_part
    A_alpha <<- sum(term1 - term2)
  }
  BC_alpha_update <- function(alpha){
    left_part <- 1/(1+K_alpha)*(z_alpha - MU) + K_alpha/((1+K_alpha)*w_alpha_star^2)*(y_alpha - MU)
    right_part <- D*partial_u_alpha %*% t(H)
    sr <- t(solve(sigma, t(right_part)))
    s <- rowSums(left_part * sr)
    B_alpha <<- sum(s)
    right_part <- (D*u_alpha - 1) %*% t(H)
    sr <- t(solve(sigma, t(right_part)))
    s <- rowSums(left_part * sr)
    C_alpha <<- sum(s)
  }




  D <- ncol(X)
  N <- nrow(X)
  H <- helmert(D)

  logX <- log(X)
  MU <- matrix(mu, nrow = N, ncol = length(mu), byrow = TRUE)
  count = 0
  while(count < 150){
    count = count + 1
    z_alpha_update(alpha)
    y_alpha_update(alpha)
    w_alpha_update(alpha)
    w_alpha_star_update(alpha)
    k_factor_update(alpha, mu = mu, sigma = sigma)
    u_alpha_update(alpha)
    partial_u_alpha_update(alpha)
    A_alpha_update(alpha)
    BC_alpha_update(alpha)
    alpha <- (B_alpha + sqrt(B_alpha^2 - 4*A_alpha*C_alpha))/(2*A_alpha)
    print(alpha)
  }
  alpha
}


alpha_est3 <- function(X, mu, sigma, alpha = 0){
  z_alpha_update <- function(alpha){
    z_alpha <<- alpha_fold(X = X, alpha = alpha, unfold = FALSE)
  }
  y_alpha_update <- function(alpha){
    y_alpha <<- alpha_fold(X = X, alpha = alpha, unfold = TRUE)
  }
  w_alpha_update <- function(alpha){
    w_alpha <<- alpha_trans(X, alpha)
  }
  w_alpha_star_update <- function(alpha){
    j_star <- apply(u_alpha, MARGIN = 1, which.min)
    w_alpha_star <<- apply(w_alpha*alpha, MARGIN = 1, min) 
  }
  k_factor_update <- function(alpha, mu, sigma){
    A <- abs(w_alpha_star)^(-2*(D-1))
    Q_1 <- 0.5*mahalanobis(z_alpha, center = mu, cov = sigma)
    Q_2 <- 0.5*mahalanobis(y_alpha, center = mu, cov = sigma)
    K_alpha <<- A*exp(Q_1 - Q_2)
  }
  u_alpha_update <- function(alpha){
    A <- alpha * logX              
    m <- apply(A, 1L, max)
    shifted_exp <- rowSums(exp(A - m))
    denom_log <- m + log(shifted_exp) 

    U <- sweep(A, 1L, denom_log, FUN = "-") 
    u_alpha <<- exp(U)
  }
  partial_u_alpha_update <- function(alpha){
    partial_u_alpha <<- (logX - rowSums(u_alpha*logX)) * u_alpha
  }
  A_alpha_update <- function(alpha){

    sr <- t(solve(sigma, t(y_alpha)))
    s <- rowSums((y_alpha - MU) * sr)

    j_star <- apply(u_alpha, MARGIN = 1, which.min)
    u_alpha_star <- apply(u_alpha, MARGIN = 1, min)
    w_part_up <- 2*D*partial_u_alpha[cbind(seq_len(N), j_star)]
    w_part_down <- D*u_alpha_star - 1
    w_part <- w_part_up/w_part_down

    term1 <- rowSums((1-D*u_alpha)*logX) 
    term2 <- K_alpha/(1+K_alpha) * (D-1 - s) * w_part
    A_alpha <<- sum(term1 - term2)
  }
  BC_alpha_update <- function(alpha){
    left_part <- 1/(1+K_alpha)*(z_alpha - MU) + K_alpha/((1+K_alpha)*w_alpha_star^2)*(y_alpha - MU)
    right_part <- D*partial_u_alpha %*% t(H)
    sr <- t(solve(sigma, t(right_part)))
    s <- rowSums(left_part * sr)
    B_alpha <<- sum(s)
    right_part <- (D*u_alpha - 1) %*% t(H)
    sr <- t(solve(sigma, t(right_part)))
    s <- rowSums(left_part * sr)
    C_alpha <<- sum(s)
  }

  D <- ncol(X)
  N <- nrow(X)
  H <- helmert(D)

  logX <- log(X)
  MU <- matrix(mu, nrow = N, ncol = length(mu), byrow = TRUE)
  count = 0
  while(count < 150){
    print(alpha)
    count = count + 1
    z_alpha_update(alpha)
    y_alpha_update(alpha)
    w_alpha_update(alpha)
    w_alpha_star_update(alpha)
    k_factor_update(alpha, mu = mu, sigma = sigma)
    u_alpha_update(alpha)
    partial_u_alpha_update(alpha)
    A_alpha_update(alpha)
    BC_alpha_update(alpha)
    partial_l <- A_alpha - B_alpha/alpha + C_alpha/alpha^2
    l <- sum(dafmvn(X, alpha = alpha, mean = mu, sigma = sigma))
    
    alpha <- alpha - l/partial_l

  }
  alpha
}

partial_l <- function(alpha, X, mu, sigma){

  k_factor_update <- function(alpha, mu, sigma){
    A <- abs(w_alpha_star)^(-2*(D-1))
    Q_1 <- 0.5*mahalanobis(z_alpha, center = mu, cov = sigma)
    Q_2 <- 0.5*mahalanobis(y_alpha, center = mu, cov = sigma)
    K_alpha <<- A*exp(Q_1 - Q_2)
  }

  partial_u_alpha_update <- function(alpha){
    partial_u_alpha <<- (logX - rowSums(u_alpha*logX)) * u_alpha
  }
  A_alpha_update <- function(alpha){

    sr <- t(solve(sigma, t(y_alpha)))
    s <- rowSums((y_alpha - MU) * sr)

    j_star <- apply(u_alpha, MARGIN = 1, which.min)
    u_alpha_star <- apply(u_alpha, MARGIN = 1, min)
    w_part_up <- 2*D*partial_u_alpha[cbind(seq_len(N), j_star)]
    w_part_down <- D*u_alpha_star - 1
    w_part <- w_part_up/w_part_down

    term1 <- rowSums((1-D*u_alpha)*logX) 
    term2 <- K_alpha/(1+K_alpha) * (D-1 - s) * w_part
    A_alpha <<- sum(term1 - term2)
  }
  BC_alpha_update <- function(alpha){
    left_part <- 1/(1+K_alpha)*(z_alpha - MU) + K_alpha/((1+K_alpha)*w_alpha_star^2)*(y_alpha - MU)
    right_part <- D*partial_u_alpha %*% t(H)
    sr <- t(solve(sigma, t(right_part)))
    s <- rowSums(left_part * sr)
    B_alpha <<- sum(s)
    right_part <- (D*u_alpha - 1) %*% t(H)
    sr <- t(solve(sigma, t(right_part)))
    s <- rowSums(left_part * sr)
    C_alpha <<- sum(s)
  }

  D <- ncol(X)
  N <- nrow(X)
  H <- helmert(D)

  logX <- log(X)
  MU <- matrix(mu, nrow = N, ncol = length(mu), byrow = TRUE)
  z_alpha <- alpha_fold(X = X, alpha = alpha, unfold = FALSE)
  y_alpha <<- alpha_fold(X = X, alpha = alpha, unfold = TRUE)
  w_alpha <<- alpha_trans(X, alpha)

  j_star <- apply(u_alpha, MARGIN = 1, which.min)
  w_alpha_star <<- apply(w_alpha*alpha, MARGIN = 1, min) 
}
