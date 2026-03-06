#' MLE of alpha-folding multivariate normal distribution given alpha
#'
#' @param X A n-by-D matrix.
#' @param alpha A numeric parameter of alpha-folding multivariate normal distribution.
#' @param tol A threshold of alpha. When |alpha| < tol, the clt is adopted.
#' @param max.iter The number of alternating iterations of solving for mu and sigma.
#' @param initial A list with field `mode` in {"default","random","user"}.
#'   - mode="default": use two deterministic starts from z_alpha_1 and z_alpha_2 moments
#'   - mode="random": use only random starts independent of data
#'       fields: n_start, mu_range, jitter, seed (optional)
#'   - mode="user": use only user-specified starts
#'       field: starts = list(list(mu=..., sigma=...), ...)
#'
#' @returns A list of estimated mu and sigma.
#' @export
afmvn_mle <- function(X,
                      alpha,
                      tol = 1e-6,
                      max.iter = 150,
                      initial = list(mode = 'default')){

  # ---------- CLT / alpha ~ 0 branch ----------
  if(abs(alpha) < tol){
    D <- ncol(X)
    H <- helmert(D)
    W <- alpha_trans(X, alpha)
    Y0 <- W %*% t(H)

    mu_hat <- colMeans(Y0)
    Sigma_hat <- cov(Y0)

    return(list(Sigma_hat = Sigma_hat, mu_hat = mu_hat))
  }

  # ---------- Precompute transformed quantities ----------
  D <- ncol(X)
  N <- nrow(X)

  H <- helmert(D)
  W <- alpha_trans(X, alpha)

  W_alpha_star <- apply(W * alpha, 1, min)
  z_alpha_1 <- W %*% t(H)
  z_alpha_2 <- z_alpha_1 / W_alpha_star^2

  # ---------- K factor ----------
  k_factor <- function(mu, sigma){
    A <- abs(W_alpha_star)^(-2 * (D - 1))
    Q_1 <- 0.5 * mahalanobis(z_alpha_1, center = mu, cov = sigma)
    Q_2 <- 0.5 * mahalanobis(z_alpha_2, center = mu, cov = sigma)
    A * exp(Q_1 - Q_2)
  }

  # ---------- Build starts according to initial$mode (mutually exclusive) ----------
  mode <- tolower(initial$mode)
  starts <- list()

  if(mode == "default"){

    mu_hat1 <- colMeans(z_alpha_1)
    Sigma_hat1 <- cov(z_alpha_1)

    mu_hat2 <- colMeans(z_alpha_2)
    Sigma_hat2 <- cov(z_alpha_2)

    starts <- list(
      list(mu = mu_hat1, sigma = Sigma_hat1),
      list(mu = mu_hat2, sigma = Sigma_hat2)
    )

  } else if(mode == "random"){

    n_start  <- if (!is.null(initial$n_start)) initial$n_start else 10
    mu_range <- if (!is.null(initial$mu_range)) initial$mu_range else 5
    jitter   <- if (!is.null(initial$jitter)) initial$jitter else 1e-6
    if (!is.null(initial$seed)) set.seed(initial$seed)

    p <- D - 1
    for (s in seq_len(n_start)) {
      mu0 <- runif(p, -mu_range, mu_range)
      A <- matrix(rnorm(p * p), p, p)
      Sigma0 <- crossprod(A) + diag(jitter, p)  # SPD
      starts[[length(starts) + 1]] <- list(mu = mu0, sigma = Sigma0)
    }

  } else if(mode == "user"){

    for(st in initial$starts){
      starts[[length(starts) + 1]] <- list(mu = st$mu, sigma = st$sigma)
    }
  } else if(mode == 'default2'){
    mu_hat1 <- colMeans(z_alpha_1)
    Sigma_hat1 <- cov(z_alpha_1)

    mu_hat2 <- colMeans(z_alpha_2)
    Sigma_hat2 <- cov(z_alpha_2)
    starts <- list(
      list(mu = (mu_hat1+mu_hat2)/2, sigma = (Sigma_hat1+Sigma_hat2)/2)
    )
  } else {
    stop("initial$mode must be one of: 'default', 'random', 'user'.")
  }

  # ---------- Fixed-point / EM-type iteration for one start ----------
  run_one_start <- function(mu_init, Sigma_init){
    mu_hat <- mu_init
    Sigma_hat <- Sigma_init

    count <- 0
    while(count <= max.iter){
      K <- k_factor(mu = mu_hat, sigma = Sigma_hat)

      w1 <- 1 / (1 + K)
      w2 <- 1 - w1

      mu_hat_updated <- colMeans(w1 * z_alpha_1 + w2 * z_alpha_2)

      z_alpha_1_shifted_weighted <- sweep(z_alpha_1, 2, mu_hat_updated, "-") * sqrt(w1)
      z_alpha_2_shifted_weighted <- sweep(z_alpha_2, 2, mu_hat_updated, "-") * sqrt(w2)

      Sigma_hat_updated <- (crossprod(z_alpha_1_shifted_weighted) +
                              crossprod(z_alpha_2_shifted_weighted)) / N

      mu_hat <- mu_hat_updated
      Sigma_hat <- Sigma_hat_updated
      count <- count + 1
    }

    list(mu_hat = mu_hat, Sigma_hat = Sigma_hat)
  }

  # ---------- Run all starts; select best by log-likelihood ----------
  best_ll <- -Inf
  best_fit <- NULL

  for(st in starts){
    fit <- run_one_start(st$mu, st$sigma)
    ll <- sum(dafmvn(X = X, alpha = alpha, mean = fit$mu_hat, sigma = fit$Sigma_hat))

    if(ll >= best_ll){
      best_ll <- ll
      best_fit <- list(mu_hat = fit$mu_hat, Sigma_hat = fit$Sigma_hat, loglik = ll)
    }
  }

  best_fit
}