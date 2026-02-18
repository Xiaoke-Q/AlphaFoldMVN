#' MCMC for AFMVN model
#'
#' Metropolis-within-Gibbs sampler for \code{dafmvn} likelihood with NIW prior.
#'
#' @param X numeric matrix (n x d)
#' @param prior list with m0, kappa0, Psi, nu
#' @param init list with alpha, mu, Sigma
#' @param warm_up,warm_band,n_iter,alpha_repropose,alpha_step,mu_step,Sigma_step tuning parameters
#' @return list with chains and acceptance rates
#' @export
afmvn_mcmc <- function(X, prior, init,
                 warm_up = 5000, warm_band = 100,
                 n_iter = 5000,
                 alpha_repropose = 10,
                 alpha_step = 0.1, mu_step = 0.1, Sigma_step = 0.1) {
  state <- list(
    alpha = init$alpha,
    mu = init$mu,
    Sigma = init$Sigma,
    loglik = loglik_af(X, init$alpha, init$mu, init$Sigma)
  )

  # warmup phase
  adaptor <- make_adaptor(
    warm_band = warm_band,
    alpha_repropose = alpha_repropose,
    init_alpha_m = 0,
    init_alpha_s = alpha_step,
    init_mu_step = mu_step,
    init_Sigma_step = Sigma_step
  )

  tuned <- adaptor$get_tuned()

  for (t in 1:warm_up) {
    out <- step_kernel(X, prior, state, tuned, alpha_repropose = alpha_repropose)
    state <- out$state
    adaptor$update(state$alpha, out$acc$alpha, out$acc$mu, out$acc$Sigma)
    tuned <- adaptor$get_tuned()
  }
  tuned_final <- tuned

  # mcmc phase
  d <- length(state$mu)
  alpha_chain <- numeric(n_iter)
  mu_chain <- matrix(NA, n_iter, d)
  Sigma_chain <- array(NA, dim = c(d, d, n_iter))

  acc_a <- 0
  acc_m <- 0
  acc_S <- 0

  for (t in 1:n_iter) {
    out <- step_kernel(X, prior, state, tuned_final, alpha_repropose = alpha_repropose)
    state <- out$state
    acc_a <- acc_a + out$acc$alpha
    acc_m <- acc_m + out$acc$mu
    acc_S <- acc_S + out$acc$Sigma

    alpha_chain[t] <- state$alpha
    mu_chain[t, ] <- state$mu
    Sigma_chain[, , t] <- state$Sigma
  }

  list(
    alpha = alpha_chain,
    mu = mu_chain,
    Sigma = Sigma_chain,
    tuned = tuned_final,
    acceptance = c(alpha = acc_a / (n_iter * alpha_repropose),
                   mu = acc_m / n_iter,
                   Sigma = acc_S / n_iter)
  )
}