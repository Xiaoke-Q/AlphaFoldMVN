step_kernel <- function(X, prior, 
                        state, tuned, 
                        alpha_repropose = 10){
  alpha <- state$alpha
  mu <- state$mu
  Sigma <- state$Sigma
  loglik_now <- state$loglik

  alpha_m <- tuned$alpha_m
  alpha_s <- tuned$alpha_s
  mu_step <- tuned$mu_step
  Sigma_step <- tuned$Sigma_step

  # alpha update
  acc_a <- 0
  for(j in 1:alpha_repropose){
    a_prop <- propose_alpha(alpha_m, alpha_s)
    loglik_prop <- loglik_af(X, a_prop, mu, Sigma)
    logr1 <- (loglik_prop - loglik_now) + r1_prior_diff(a_prop, alpha, alpha_m, alpha_s)

    if(is.finite(logr1) && (log(runif(1)) < logr1)){
      alpha <- a_prop
      loglik_now <- loglik_prop
      acc_a <- acc_a + 1
    }
  }

  #mu updtae
  acc_m <- 0
  mu_prop <- propose_mu(mu, mu_step, Sigma)
  loglik_prop <- loglik_af(X, alpha, mu_prop, Sigma)
  logr2 <- (loglik_prop - loglik_now) + r2_prior_diff(mu_prop, mu, Sigma, prior)

  if(is.finite(logr2) && (log(runif(1)) < logr2)){
    mu <- mu_prop
    loglik_now <- loglik_prop
    acc_m <- 1
  }

  # sigma update
  acc_S <- 0
  Sigma_prop <- propose_Sigma(Sigma, Sigma_step)
  loglik_prop <- loglik_af(X, alpha, mu, Sigma_prop)
  logr3 <- (loglik_prop - loglik_now) + r3_prior_diff(Sigma_prop, Sigma, mu, prior)

  if(is.finite(logr3) && (log(runif(1)) < logr3)){
    Sigma <- Sigma_prop
    loglik_now <- loglik_prop
    acc_S <- 1
  }

  list(
    state = list(alpha = alpha, mu = mu, Sigma = Sigma, loglik = loglik_now),
    acc = list(alpha = acc_a, mu = acc_m, Sigma = acc_S)
  )
}