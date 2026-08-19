devtools::load_all()

set.seed(1)

alpha_true <- 0.7
eta_true <- atanh(alpha_true)

mean_true <- c(-1, 1 ,0, 0.5, 0.7)
sigma_true <- diag(5)

sim <- rafmvn(
  n = 200,
  alpha = alpha_true,
  mean = mean_true,
  sigma = sigma_true
)

X <- sim$X
folded <- sim$folded

N <- nrow(X)
p <- ncol(X) - 1

m_eta_init <- 0.3
s2_eta_init <- 0.1

m_q_init <- rep(0, p)

nu_q_init <- p + 20
Psi_q_init <- diag(p) / nu_q_init

m_0 <- rep(0, p)
kappa_0 <- 0.01
Psi_0 <- diag(p)
nu_0 <- p + 2

fit <- afmvn_vba(
  X = X,
  m_eta_init = m_eta_init,
  s2_eta_init = s2_eta_init,
  m_q_init = m_q_init,
  Psi_q_init = Psi_q_init,
  nu_q_init = nu_q_init,
  m_0 = m_0,
  kappa_0 = kappa_0,
  Psi_0 = Psi_0,
  nu_0 = nu_0,
  a_0 = 0,
  s0_sq = 10,
  n_quad = 15,
  eps_alpha = 1e-8,
  eta_control = list(maxit = 150L),
  max_iter = 200, 
  tol = 1e-5
)

fit$trace
head(fit$trace)

alpha_true
tanh(fit$m_eta)
eta_true
fit$m_eta
fit$s2_eta





plot(fit$trace[,1],fit$trace[,3])
abline(h=alpha_true, lty=2)

plot(fit$trace[,1], log(unlist(fit$trace[,4])), ylab = "log(s2_eta)")

plot(fit$trace[,1], fit$trace[,9], ylab = "max_change")

plot(fit$trace[,1], fit$trace[,5], ylab = "eta_objective")





set.seed(1027)
alpha_values <- c(-0.7, -0.3, 0,0.3, 0.7)

n_rep <- 20
N <- 200
p <- 5

alpha_hat <- matrix(NA, n_rep, length(alpha_values))
mean_hat <- list()
sigma_hat <- list()

mean_values <- list()
sigma_values <- list()
fit_values <- list()

for (b in 1:n_rep) {
  Lambda_true <- rWishart(1, p + 10, diag(p)/(9))[,,1]
  sigma_true <- solve(Lambda_true)
  mean_true <- rnorm(p)

  mean_values[[b]] <- mean_true
  sigma_values[[b]] <- sigma_true
  mean_hat[[b]] <- list()
  sigma_hat[[b]] <- list()
  fit_values[[b]] <- list()

  for (j in 1:length(alpha_values)) {
    alpha_true <- alpha_values[j]

    sim <- rafmvn(n = N, alpha = alpha_true, mean = mean_true, sigma = sigma_true)
    X <- sim$X

    m_eta_init <- 0
    s2_eta_init <- 0.1
    Z_init <- alpha_fold(X = X, alpha = tanh(m_eta_init), unfold = FALSE)
    m_q_init <- colMeans(Z_init)
    S_init <- cov(Z_init)
    nu_q_init <- p + 20
    Psi_q_init <- solve(S_init) / nu_q_init
    cat("b = ", b, ", j = ", j, "\n")
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
    fit_values[[b]][[j]] <- fit
  }
}

save.image('dev/sim1.RData')
load('dev/sim1.RData')
alpha_hat <- matrix(NA, n_rep, length(alpha_values))
mean_hat <- list()
sigma_hat <- list()

for (b in 1:n_rep) {
  mean_hat[[b]] <- list()
  sigma_hat[[b]] <- list()

  for (j in 1:length(alpha_values)) {
    fit <- fit_values[[b]][[j]]

    alpha_hat[b, j] <- tanh(fit$m_eta)
    mean_hat[[b]][[j]] <- fit$m_q
    sigma_hat[[b]][[j]] <- solve(fit$Psi_q)/(fit$nu_q - p - 1)
  }
}

alpha_error <- alpha_hat - rep(alpha_values, each = n_rep)
mean_error <- lapply(1:n_rep, function(b) lapply(1:length(alpha_values), function(j) mean_hat[[b]][[j]] - mean_values[[b]]))
sigma_error <- lapply(1:n_rep, function(b) lapply(1:length(alpha_values), function(j) sigma_hat[[b]][[j]] - sigma_values[[b]]))
alpha_error
mean_error
sigma_error

alpha_bias <- colMeans(alpha_error)
alpha_rmse <- sqrt(colMeans(alpha_error^2))
alpha_rmse

mean_bias <- lapply(1:length(alpha_values), function(j) Reduce("+", lapply(mean_error, function(x) x[[j]]))/n_rep)
mean_rmse <- lapply(1:length(alpha_values), function(j) sqrt(Reduce("+", lapply(mean_error, function(x) x[[j]]^2))/n_rep))
mean_rmse

sigma_bias <- lapply(1:length(alpha_values), function(j) Reduce("+", lapply(sigma_error, function(x) x[[j]]))/n_rep)
sigma_rmse <- lapply(1:length(alpha_values), function(j) sqrt(Reduce("+", lapply(sigma_error, function(x) x[[j]]^2))/n_rep))
sigma_rmse
