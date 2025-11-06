unfold <- function(W, alpha){
  W_alpha_star <- apply(alpha*W, MARGIN = 1, min)
  W <- W/(W_alpha_star)^2
  W
}

unfold_2 <- function(W, alpha){
  W_alpha_star <- alpha*apply(W, MARGIN = 1, min)
  W <- W/(W_alpha_star)^2
  W
}