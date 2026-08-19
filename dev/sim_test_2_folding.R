devtools::load_all()

set.seed(1)

alpha_values <- c(0.7)
N_values <- c(100, 200, 500)
p_values <- c(2, 10, 20)
n_rep <- 10

mean_values <- list()
sigma_values <- list()
fit_values <- list()

for (d in 1:length(p_values)) {
  p <- p_values[d]

  mean_values[[d]] <- list()
  sigma_values[[d]] <- list()
  fit_values[[d]] <- list()

  for (b in 1:n_rep) {
    Lambda_true <- rWishart(1, p + 10, diag(p)/9)[,,1]
    sigma_true <- solve(Lambda_true)
    mean_true <- rnorm(p)

    mean_values[[d]][[b]] <- mean_true
    sigma_values[[d]][[b]] <- sigma_true
  }

  for (s in 1:length(N_values)) {
    N <- N_values[s]
    fit_values[[d]][[s]] <- list()

    for (b in 1:n_rep) {
      mean_true <- mean_values[[d]][[b]]
      sigma_true <- sigma_values[[d]][[b]]
      fit_values[[d]][[s]][[b]] <- list()

      for (j in 1:length(alpha_values)) {
        alpha_true <- alpha_values[j]

        sim <- rafmvn(n = N, alpha = alpha_true, mean = mean_true, sigma = sigma_true)
        X <- sim$X

        m_eta_init <- 0.3
        s2_eta_init <- 0.1
        Z_init <- alpha_fold(X = X, alpha = tanh(m_eta_init), unfold = FALSE)
        m_q_init <- colMeans(Z_init)
        S_init <- cov(Z_init)
        nu_q_init <- p + 20
        Psi_q_init <- solve(S_init)/nu_q_init

        fit <- afmvn_vba(
          X = X,
          m_eta_init = m_eta_init,
          s2_eta_init = s2_eta_init,
          m_q_init = m_q_init,
          Psi_q_init = Psi_q_init,
          nu_q_init = nu_q_init,
          m_0 = rep(0, p),
          kappa_0 = 0.01,
          Psi_0 = diag(p),
          nu_0 = p + 2,
          a_0 = 0,
          s0_sq = 10,
          n_quad = 15,
          eps_alpha = 1e-8,
          eta_control = list(maxit = 150),
          max_iter = 200, 
          tol = 1e-5
        )

        fit_values[[d]][[s]][[b]][[j]] <- fit
      }
    }
  }
}

save.image('dev/sim2.RData')
load('dev/sim2.RData')

alpha_hat <- list()
mean_hat <- list()
sigma_hat <- list()

for (d in 1:length(p_values)) {
  p <- p_values[d]

  alpha_hat[[d]] <- list()
  mean_hat[[d]] <- list()
  sigma_hat[[d]] <- list()

  for (s in 1:length(N_values)) {
    alpha_hat[[d]][[s]] <- matrix(NA, n_rep, length(alpha_values))
    mean_hat[[d]][[s]] <- list()
    sigma_hat[[d]][[s]] <- list()

    for (b in 1:n_rep) {
      mean_hat[[d]][[s]][[b]] <- list()
      sigma_hat[[d]][[s]][[b]] <- list()

      for (j in 1:length(alpha_values)) {
        fit <- fit_values[[d]][[s]][[b]][[j]]

        alpha_hat[[d]][[s]][b, j] <- tanh(fit$m_eta)
        mean_hat[[d]][[s]][[b]][[j]] <- fit$m_q
        sigma_hat[[d]][[s]][[b]][[j]] <- solve(fit$Psi_q)/(fit$nu_q - p - 1)
      }
    }
  }
}

d <- which(p_values == 20)
s <- which(N_values == 200)
b <- 1
j <- which(alpha_values == 0.7)

alpha_hat[[d]][[s]][b, j]
mean_values[[d]][[b]]
mean_hat[[d]][[s]][[b]][[j]]
sigma_values[[d]][[b]]
sigma_hat[[d]][[s]][[b]][[j]]


alpha_error <- list()
mean_error <- list()
sigma_error <- list()

for (d in 1:length(p_values)) {
  alpha_error[[d]] <- list()
  mean_error[[d]] <- list()
  sigma_error[[d]] <- list()

  for (s in 1:length(N_values)) {
    alpha_error[[d]][[s]] <- alpha_hat[[d]][[s]] - rep(alpha_values, each = n_rep)
    mean_error[[d]][[s]] <- list()
    sigma_error[[d]][[s]] <- list()

    for (b in 1:n_rep) {
      mean_error[[d]][[s]][[b]] <- list()
      sigma_error[[d]][[s]][[b]] <- list()

      for (j in 1:length(alpha_values)) {
        mean_error[[d]][[s]][[b]][[j]] <- mean_hat[[d]][[s]][[b]][[j]] - mean_values[[d]][[b]]
        sigma_error[[d]][[s]][[b]][[j]] <- sigma_hat[[d]][[s]][[b]][[j]] - sigma_values[[d]][[b]]
      }
    }
  }
}
alpha_rmse <- list()
mean_rmse <- list()
sigma_rmse <- list()

for (d in 1:length(p_values)) {
  alpha_rmse[[d]] <- list()
  mean_rmse[[d]] <- list()
  sigma_rmse[[d]] <- list()

  for (s in 1:length(N_values)) {
    alpha_rmse[[d]][[s]] <- sqrt(colMeans(alpha_error[[d]][[s]]^2))

    mean_rmse[[d]][[s]] <- lapply(1:length(alpha_values), function(j) {
      sqrt(Reduce("+", lapply(mean_error[[d]][[s]], function(x) x[[j]]^2))/n_rep)
    })

    sigma_rmse[[d]][[s]] <- lapply(1:length(alpha_values), function(j) {
      sqrt(Reduce("+", lapply(sigma_error[[d]][[s]], function(x) x[[j]]^2))/n_rep)
    })
  }
}


d <- which(p_values == 2)
s <- which(N_values == 200)

alpha_rmse[[d]][[s]]
mean_rmse[[d]][[s]]
sigma_rmse[[d]][[s]]

tab <- do.call(rbind, 
  lapply(1:length(p_values), 
  function(d) do.call(rbind, 
    lapply(1:length(N_values), function(s) alpha_rmse[[d]][[s]]))))

rownames(tab) <- as.vector(sapply(p_values, function(p) paste0("p=", p, ", N=", N_values)))
colnames(tab) <- paste0("alpha=", alpha_values)

round(tab, 4)


d <- which(p_values == 20)
j <- which(alpha_values == 0.7)

cbind(
  N = N_values,
  mean = sapply(1:length(N_values), function(s) mean(alpha_hat[[d]][[s]][,j])),
  bias = sapply(1:length(N_values), function(s) mean(alpha_hat[[d]][[s]][,j]) - alpha_values[j]),
  sd = sapply(1:length(N_values), function(s) sd(alpha_hat[[d]][[s]][,j])),
  rmse = sapply(1:length(N_values), function(s) sqrt(mean((alpha_hat[[d]][[s]][,j] - alpha_values[j])^2)))
)

s <- which(N_values == 500)
b <- which.max(abs(alpha_hat[[d]][[s]][,j] - alpha_values[j]))

alpha_hat[[d]][[s]][b,j]
fit_values[[d]][[s]][[b]][[j]]$trace
