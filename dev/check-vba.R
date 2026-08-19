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
  eta_control = list(maxit = 150L)
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


m_grid <- seq(fit$m_eta - 1.5, fit$m_eta + 1.5, length.out = 60L)
log_s_grid <- seq(log(0.005), log(1), length.out = 60L)

surf <- eta_objective_surface(
  X = X, r = fit$r,
  m_q = fit$m_q, Psi_q = fit$Psi_q, nu_q = fit$nu_q,
  m_grid = m_grid, log_s_grid = log_s_grid,
  a0 = 0, s0_sq = 10,
  n_quad = 20, eps_alpha = 1e-8
)

plot(surf_final, nlevels = 60)
abline(v = eta_true, lty = 2)




alpha_starts <- c(-0.7, -0.3, 0, 0.3, 0.7)

fits <- vector("list", 5)

for (k in seq_along(alpha_starts)) {
  alpha0 <- alpha_starts[k]

  comp <- .afmvn_components_alpha(
    X = X,
    alpha = alpha0,
    eps_alpha = 1e-8
  )

  Z <- comp$z0
  S <- stats::cov(Z)

  fits[[k]] <- afmvn_vba(
    X = X,
    m_eta_init = atanh(alpha0),
    s2_eta_init = 0.25,
    m_q_init = colMeans(Z),
    Psi_q_init = solve(S) / (p + 20),
    nu_q_init = p + 20,
    m_0 = rep(0, p),
    kappa_0 = 0.01,
    Psi_0 = diag(p),
    nu_0 = p + 2
  )
}

sapply(fits, function(fit) tail(fit$trace$max_change, 1))
sapply(fits, function(fit) tanh(fit$m_eta))



starts <- expand.grid(
  alpha0 = c(-0.7, -0.3, 0, 0.3, 0.7),
  s2_eta0 = c(0.05, 0.25, 0.75)
)

nu_q_init <- p + 20
Psi_q_init <- diag(p)/nu_q_init

m_0 <- rep(0, p)
kappa_0 <- 0.01
Psi_0 <- diag(p)
nu_0 <- p + 2

fits <- vector("list", nrow(starts))
.afmvn_init_from_alpha <- function(X, alpha0,
                                   s2_eta_init = 0.25,
                                   nu_extra = 20L,
                                   ridge = 1e-6,
                                   eps_alpha = 1e-8) {
  X <- as.matrix(X)
  p <- ncol(X) - 1

  comp <- .afmvn_components_alpha(
    X = X,
    alpha = alpha0,
    eps_alpha = eps_alpha
  )

  Z <- comp$z0

  m_q_init <- colMeans(Z)

  S <- cov(Z)
  S <- 0.5 * (S + t(S))
  S <- S + ridge * diag(p)

  nu_q_init <- p + nu_extra
  Psi_q_init <- solve(S) / nu_q_init
  Psi_q_init <- 0.5*(Psi_q_init + t(Psi_q_init))

  list(
    alpha = alpha0,
    m_eta = atanh(alpha0),
    s2_eta = s2_eta_init,
    m_q = m_q_init,
    Psi_q = Psi_q_init,
    nu_q = nu_q_init,
    Z = Z
  )
}


for (k in seq_len(nrow(starts))) {
  alpha0 <- starts$alpha0[k]
  s2_eta0 <- starts$s2_eta0[k]

  init <- .afmvn_init_from_alpha(
    X = X,
    alpha0 = alpha0,
    s2_eta_init = s2_eta0,
    nu_extra = 20L
  )

  fits[[k]] <- afmvn_vba(
    X = X,
    m_eta_init = init$m_eta,
    s2_eta_init = init$s2_eta,
    m_q_init = init$m_q,
    Psi_q_init = init$Psi_q,
    nu_q_init = init$nu_q,
    m_0 = m_0,
    kappa_0 = kappa_0,
    Psi_0 = Psi_0,
    nu_0 = nu_0
  )
}


alpha_hat <- sapply(fits, function(f) tanh(f$m_eta))
alpha0
alpha_hat
mean(alpha_hat)
sd(alpha_hat)


s2_eta_hat <- sapply(fits, function(f) f$s2_eta)
s2_eta_hat


mu_hat <- lapply(fits, function(f) f$m_q)
mean_true
mu_hat

sigma_hat <- lapply(fits, function(f) {
  solve(f$nu_q * f$Psi_q)
})
sigma_true
sigma_hat
