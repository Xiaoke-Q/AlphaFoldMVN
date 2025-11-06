#' Alpha transformation of compositional data
#'
#' @param X A numeric vector, matrix, or data frame of positive values. 
#' If a vector is supplied, it is treated as a single row.
#' If matrix or data frame are supplied, they are row stochastic.
#' @param alpha A numeric scalar specifying the transformation parameter.
#' Values close to 0 (<1e-8) adpot the log-ratio transform.
#' @param unfold Logical; Transform to outside of the simplex when unfold = TRUE.
#' @param check_one Logical; Should we check X is row stochastic?
#' 
#' @param renorm Logical: Should we normalize rows of X?
#' @returns A matrix with the same dimention of X.
#'
#' @export
#' @examples X <- c(0.1, 0.6, 0.3) 
#' alpha_trans(X = X , alpha = 0.7)
alpha_trans <- function(X, alpha, unfold = FALSE, check_one = TRUE, renorm = FALSE){
  is_vec <- is.null(dim(X))
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is_vec) X <- matrix(X, nrow = 1L)
  
  # check positivity
  if (any(!(X > 0))) {
    stop("alpha_trans: X must be strictly positive.")
  }
  # check row stochastic
  rs <- rowSums(X)
  if (check_one) {
    bad <- abs(rs - 1) > 1e-8
    if (any(bad)) {
      if (!renorm) {
        stop("alpha_trans: rows must sum to 1 (row-stochastic). ",
             "Found rows with sums not equal to 1. ",
             "Set renorm=TRUE to rescale rows to the simplex.")
      } else {
        X <- X/rs
      }
    }
  } else if (renorm) {
    X <- X/rs
  }
  
  logX <- log(X)

  if (abs(alpha) <= 1e-8) {
    mu <- rowMeans(logX)
    W <- sweep(logX, 1L, mu, FUN = "-")
  }else{
    A <- alpha * logX              
    m <- apply(A, 1L, max)
    shifted_exp <- rowMeans(exp(A - m))
    denom_log <- m + log(shifted_exp) 

    U <- sweep(A, 1L, denom_log, FUN = "-") 
    W <- (exp(U) - 1) / alpha
    if(length(unfold) == 1) unfold = rep(unfold, nrow(X))
    W_alpha_star <- apply(alpha*W[unfold,], MARGIN = 1, min)
    W[unfold,] <- W[unfold,]/(W_alpha_star)^2
  }
  W
}



#' The inverse of alpha_trans
#'
#' @param W A numeric vector or matrix. 
#' If a vector is supplied, it is treated as a single row.
#' The row sums should be 0.
#' @param alpha A numeric scalar specifying the transformation parameter.
#' Values close to 0 (<1e-8) adpot the log-ratio transform.
#'
#' @returns A numeric vector or matrix, which is the inverse of W given alpha.
#'
#' @export
#' @examples alpha_trans_inverse(alpha_trans(X = c(0.1, 0.6, 0.3), alpha = 0.7), alpha = 0.7)
alpha_trans_inverse <- function(W, alpha) {
  is_vec <- is.null(dim(W))
  if (is_vec) W <- matrix(W, nrow = 1L)
  
  ws <- rowSums(W)
  if(any(abs(ws) > 1e-8)) stop("alpha_invers: the row sums of W should be zero.")

  D <- ncol(W)
  N <- nrow(W)
  
  if (abs(alpha) <= 1e-8) {
    tmp <- exp(W)
    folded <- rep(0, nrow(W))
  } else {
    lb <- min(-1/alpha, (D-1)/alpha)
    rb <- max(-1/alpha, (D-1)/alpha)
    inside_index <- rowSums((W<lb)|(W>rb)) == 0
    q_alpha_star <- rep(NA, nrow = N)
    q_alpha_star[!inside_index] <- 
      apply(W[!inside_index, ,drop = FALSE] * alpha, MARGIN = 1, min)
    q_alpha_star[inside_index] <- 1
    M <- W/q_alpha_star^2

    base <- 1 + alpha * M
    base[base < 0 & base > -sqrt(.Machine$double.eps)] <- 0
    if (any(base < 0)) {
      stop("alpha_inverse: encountered 1 + alpha*W < 0.")
    }
    tmp <- base^(1/alpha)
    folded <- !inside_index
  }

  X <- tmp / rowSums(tmp)
  return(list(X = X, folded = folded))
}