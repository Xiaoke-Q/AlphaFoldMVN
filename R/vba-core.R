.afmvn_quad_context <- function(m_q, Psi_q) {
  m_q <- as.numeric(m_q)
  Psi_q <- as.matrix(Psi_q)

  Psi_m <- as.numeric(Psi_q %*% m_q)

  list(
    Psi_q = Psi_q,
    Psi_m = Psi_m,
    c_q = drop(crossprod(m_q, Psi_m)),
    vec_Psi_t = as.vector(t(Psi_q))
  )
}



.afmvn_branch_scores <- function(mbj, quad, nu_q) {
  mbj$J - 0.5 * nu_q * quad
}

.afmvn_row_softmax <- function(score) {
  score <- as.matrix(score)

  row_max <- apply(score, 1, max)

  exp_score <- exp(score - row_max)
  exp_score / rowSums(exp_score)
}

.afmvn_quad_rows <- function(z, quad_ctx) {
  rowSums((z %*% quad_ctx$Psi_q) * z) - 2 * drop(z %*% quad_ctx$Psi_m) + quad_ctx$c_q
}

.afmvn_update_r <- function(J, quad, nu_q) {
  score <- J - 0.5 * nu_q * quad
  r <- .afmvn_row_softmax(score)

  list(r = r, score = score)
}