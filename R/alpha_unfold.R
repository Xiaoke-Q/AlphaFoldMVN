alpha_unfold <- function(X, alpha, unfold = FALSE){
  if(alpha == 0){
    W <- log(X) - rowMeans(log(X))
  }else{
    u <- X^alpha/rowMeans(X^alpha)
    W <- (u-1)/alpha
    if(unfold){
      W_alpha_star <- apply(alpha*W, MARGIN = 1, min)
      W <- W/(W_alpha_star)^2
    }
  }
  W
}