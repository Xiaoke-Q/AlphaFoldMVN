#' Alpha folded transfomration of compositional data
#'
#' @param X A numeric vector, matrix, or data frame of positive values. 
#' If a vector is supplied, it is treated as a single row.
#' If matrix or data frame are supplied, they are row stochastic.
#' @param alpha A numeric scalar specifying the transformation parameter.
#' Values close to 0 (<1e-8) adpot the log-ratio transform.
#' @param unfold Logical; Transform to outside of the simplex when unfold = TRUE.
#' @param check_one Logical; Should we check X is row stochastic?
#' 
#'
#' @param renorm Logical: Should we normalize rows of X?
#' @returns A matrix with the same dimention of X.
#'
#' @export
#' @examples X <- c(0.1, 0.6, 0.3) 
#' alpha_fold(X = X , alpha = 0.7)
alpha_fold <- function(X, 
                       alpha, 
                       unfold = FALSE, 
                       check_one = TRUE, 
                       renorm = FALSE){
  is_vec <- is.null(dim(X))
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is_vec) X <- matrix(X, nrow = 1L)
  
  N <- nrow(X)
  D <- ncol(X)
  H <- helmert(D)
  Y <- alpha_trans(X = X, 
                   alpha = alpha, 
                   unfold = unfold, 
                   check_one = check_one, 
                   renorm = renorm)%*%t(H)
  Y
}

#' Inverse of alpha fold transformation
#'
#' @param Y A numeric vector or matrix. 
#' If a vector is supplied, it is treated as a single row.
#' @param alpha A numeric scalar specifying the transformation parameter.
#' Values close to 0 (<1e-8) adpot the log-ratio transform.
#'
#' @returns A numeric vector or matrix.
#'
#' @export
#' @examples alpha_fold_inverse(Y = c(0.1, 0.2), alpha = 0.1)
alpha_fold_inverse <- function(Y, alpha){
  is_vec <- is.null(dim(Y))
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  if (is_vec) Y <- matrix(Y, nrow = 1L)
  
  D <- ncol(Y) + 1
  H <- helmert(D)

  alpha_trans_inverse(Y%*%H, alpha = alpha)
}


