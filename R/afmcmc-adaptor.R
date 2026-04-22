make_adaptor <- function(warm_band, alpha_repropose, 
                         init_alpha_m = 0, init_alpha_s = 0.2,
                         init_mu_step = 0.1, init_Sigma_step = 0.1, 
                         min_alpha_s = 0.1){
  buf_a <- numeric(warm_band)
  buf_m <- numeric(warm_band)
  buf_S <- numeric(warm_band)
  sum_a <- 0
  sum_m <- 0
  sum_S <- 0

  # slide windows for alpha
  buf_z <- numeric(warm_band)
  sum_z <- 0
  sum_z2 <- 0
  filled <- 0

  t <- 0
  alpha_m <- init_alpha_m
  alpha_s <- max(init_alpha_s, min_alpha_s)
  mu_step <- init_mu_step
  Sigma_step <- init_Sigma_step

  list(
    get_tuned = function() list(alpha_m = alpha_m, 
                                alpha_s = alpha_s, 
                                mu_step = mu_step, 
                                Sigma_step = Sigma_step),

    update = function(alpha, acc_a, acc_m, acc_S){
      t <<- t + 1
      pos <- (t - 1)%%warm_band + 1

      sum_a <<- sum_a - buf_a[pos] + acc_a
      buf_a[pos] <<- acc_a
      sum_m <<- sum_m - buf_m[pos] + acc_m
      buf_m[pos] <<- acc_m
      sum_S <<- sum_S - buf_S[pos] + acc_S
      buf_S[pos] <<- acc_S

      den <- if(t < warm_band) t else warm_band
      rate_a <- (sum_a/den) / alpha_repropose
      rate_m <- (sum_m/den)
      rate_S <- (sum_S/den)

      z <- alpha_to_z(alpha)
      if (filled < warm_band) {
        filled <<- filled + 1
      }else{
        z_old <- buf_z[pos]
        sum_z <<- sum_z - z_old
        sum_z2 <<- sum_z2 - z_old * z_old
      }
      buf_z[pos] <<- z
      sum_z <<- sum_z + z
      sum_z2 <<- sum_z2 + z * z

      if(pos == warm_band){
        m_hat <- sum_z / filled
        v_hat <- sum_z2 / filled - m_hat * m_hat
        v_hat <- max(v_hat, 1e-6)
        alpha_m <<- as.numeric(m_hat)
        alpha_s <<- max(as.numeric(sqrt(v_hat)), min_alpha_s)

        if (rate_a < 0.15) alpha_s <<- max(alpha_s * 0.7, min_alpha_s)
        if (rate_a > 0.45) alpha_s <<- max(alpha_s * 1.3, min_alpha_s)

        if (rate_m < 0.15) mu_step <<- mu_step * 0.7
        if (rate_m > 0.45) mu_step <<- mu_step * 1.3
        if (rate_S < 0.15) Sigma_step <<- Sigma_step * 0.7
        if (rate_S > 0.45) Sigma_step <<- Sigma_step * 1.3
      }

      invisible(NULL)
    }
  )
}