devtools::load_all()

set.seed(1)
alpha_true <- 0.7
eta_true <- atanh(alpha_true)
mean_true <- c(-1, 1)
sigma_true <- diag(2)

sim <- rafmvn(n = 200, alpha = alpha_true, 
  mean = mean_true,
  sigma = sigma_true)

X <- sim$X
folded <- sim$folded

p <- ncol(X) - 1
m_eta <- eta_true
s2_eta <- 0.05

m_q <- mean_true

nu_q <- p + 200
Psi_q <- solve(sigma_true) / nu_q


r <- cbind(as.numeric(!folded), as.numeric(folded))
r <- cbind(
  ifelse(folded, 0.20, 0.80),
  ifelse(folded, 0.80, 0.20)
)

mbj <- .afmvn_components_qeta(
  X = X,
  m_eta = m_eta,
  s2_eta = s2_eta,
  n_quad = 20,
  eps_alpha = 1e-8
)

out_r <- .afmvn_update_r(
  M = mbj$M,
  B = mbj$B,
  J = mbj$J,
  m_q = m_q,
  Psi_q = Psi_q,
  nu_q = nu_q
)

r <- out_r$r


m_grid <- seq(eta_true - 2, eta_true + 2, length.out = 80)
log_s_grid <- seq(log(0.005), log(2), length.out = 80)

surf <- eta_objective_surface(
  X = X, r = r,
  m_q = m_q, Psi_q = Psi_q, nu_q = nu_q,
  m_grid = m_grid, log_s_grid = log_s_grid,
  a0 = 0, s0_sq = 10,
  n_quad = 20, eps_alpha = 1e-8
)

plot(surf, nlevels = 100, relative = FALSE)

abline(v = eta_true, lty = 2)
abline(h = 0.5*log(s2_eta), lty = 2)
points(eta_true, 0.5*log(s2_eta), pch = 19, col = 'green' )

obj <- function(par) {
  m_eta_cur <- par[1]
  log_s_cur <- par[2]
  s2_eta_cur <- exp(2*log_s_cur)
  res <- .afmvn_eta_objective(
    X = X, 
    m_eta = m_eta_cur, s2_eta = s2_eta_cur, 
    r = r,
    m_q = m_q, Psi_q = Psi_q, nu_q = nu_q,
    a0 = 0, s0_sq = 10,
    n_quad = 20, eps_alpha = 1e-8) 

  -res
}
opt <- optim(
  par = c(m_eta, log(sqrt(s2_eta))),
  fn = obj,
  method = "BFGS"
)
opt <- optim(
  par = c(m_eta, log(sqrt(s2_eta))),
  fn = obj,
  method = "Nelder-Mead"
)
m_eta_hat <- opt$par[1]
s2_eta_hat <- exp(2*opt$par[2])
L_hat <- -opt$value

points(m_eta_hat, 0.5*log(s2_eta_hat), pch = 19, col = 'blue')



k_bad <- which.min(surf$L[1,])
k_bad
surf$log_s_grid[k_bad]
surf$s_grid[k_bad]
surf$L[k_bad, 1]

m_bad <- surf$m_grid[1]
s2_bad <- surf$s2_grid[k_bad]

q_bad <- .eta_quadrature(
  m_eta = m_bad,
  s2_eta = s2_bad,
  n_quad = 20
)

alpha_bad <- tanh(q_bad$eta)

cbind(
  node = seq_along(q_bad$eta),
  eta = q_bad$eta,
  alpha = alpha_bad,
  weight = q_bad$weight
)



.afmvn_eta_objective_parts <- function(X, m_eta, s2_eta, r,
                                       m_q, Psi_q, nu_q,
                                       a0 = 0, s0_sq = 10,
                                       n_quad = 20, eps_alpha = 1e-8) {
  X <- as.matrix(X)
  r <- as.matrix(r)
  Psi_q <- as.matrix(Psi_q)

  N <- nrow(X)

  mbj <- .afmvn_components_qeta(
    X = X,
    m_eta = m_eta,
    s2_eta = s2_eta,
    n_quad = n_quad,
    eps_alpha = eps_alpha
  )

  M <- mbj$M
  B <- mbj$B
  J <- mbj$J

  KL_term <- .afmvn_eta_prior_KL(
    m_eta = m_eta,
    s2_eta = s2_eta,
    a0 = a0,
    s0_sq = s0_sq
  )

  j_by_branch <- colSums(r * J)
  quad_by_branch <- numeric(2)

  for (g in 1:2) {
    for (i in seq_len(N)) {
      M_ig <- M[i, g, ]
      B_ig <- B[i, g, , ]

      quad_ig <- sum(Psi_q * t(B_ig)) - 2*as.numeric(t(m_q) %*% Psi_q %*% M_ig) + as.numeric(t(m_q) %*% Psi_q %*% m_q)

      quad_by_branch[g] <- quad_by_branch[g] + r[i,g]*quad_ig
    }
  }

  quad_contribution_by_branch <- -0.5*nu_q*quad_by_branch

  list(
    value = KL_term + sum(j_by_branch) + sum(quad_contribution_by_branch),
    KL_term = KL_term,
    j_by_branch = j_by_branch,
    quad_by_branch = quad_by_branch,
    quad_contribution_by_branch = quad_contribution_by_branch,
    near_zero_node = mbj$near_zero_node,
    eta = mbj$eta,
    alpha = mbj$alpha
  )
}

k_bad <- which.min(surf$L[, 1])
k_ok <- k_bad - 1

m_val <- surf$m_grid[1]

bad <- .afmvn_eta_objective_parts(
  X = X,
  m_eta = m_val,
  s2_eta = surf$s2_grid[k_bad],
  r = r,
  m_q = m_q,
  Psi_q = Psi_q,
  nu_q = nu_q,
  n_quad = 20,
  eps_alpha = 1e-8
)

ok <- .afmvn_eta_objective_parts(
  X = X,
  m_eta = m_val,
  s2_eta = surf$s2_grid[k_ok],
  r = r,
  m_q = m_q,
  Psi_q = Psi_q,
  nu_q = nu_q,
  n_quad = 20,
  eps_alpha = 1e-8
)

ok[c("value", "KL_term", "j_by_branch", "quad_by_branch", "quad_contribution_by_branch")]
bad[c("value", "KL_term", "j_by_branch", "quad_by_branch", "quad_contribution_by_branch")]


mbj_bad <- .afmvn_components_qeta(
  X = X,
  m_eta = m_val,
  s2_eta = surf$s2_grid[k_bad],
  n_quad = 20,
  eps_alpha = 1e-8
)

quad_i1 <- numeric(nrow(X))

for (i in seq_len(nrow(X))) {
  M_ig <- mbj_bad$M[i, 2, ]
  B_ig <- mbj_bad$B[i, 2, , ]

  quad_i1[i] <- sum(Psi_q * t(B_ig)) -
    2 * as.numeric(t(m_q) %*% Psi_q %*% M_ig) +
    as.numeric(t(m_q) %*% Psi_q %*% m_q)
}

ord <- order(quad_i1, decreasing = TRUE)

head(data.frame(
  i = ord,
  r_i1 = r[ord, 2],
  quad_i1 = quad_i1[ord],
  weighted_quad_i1 = r[ord, 2] * quad_i1[ord]
), 10)


weighted_quad_i1 <- r[, 2]*quad_i1

ord_w <- order(weighted_quad_i1, decreasing = TRUE)

head(data.frame(
  i = ord_w,
  r_i1 = r[ord_w, 2],
  quad_i1 = quad_i1[ord_w],
  weighted_quad_i1 = weighted_quad_i1[ord_w]
), 20)
