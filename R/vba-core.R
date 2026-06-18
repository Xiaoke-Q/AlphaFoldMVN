.afmvn_branch_scores <- function(M, B, J, m_q, Psi_q, nu_q) {
  N <- dim(M)[1]
  p <- dim(M)[3]

  score <- matrix(NA, nrow = N, ncol = 2)
  c_q <- as.numeric(t(m_q)%*%Psi_q%*%m_q)
  v_q <- as.numeric(Psi_q%*%m_q)

  for (g in 1:2) {
    for (i in seq_len(N)) {
      M_ig <- M[i, g, ]
      B_ig <- B[i, g, , ]

      quad_ig <- sum(Psi_q * t(B_ig)) - 2*sum(v_q*M_ig) + c_q
      score[i, g] <- J[i, g] - 0.5*nu_q*quad_ig
    }
  }
  score
}

.afmvn_row_softmax <- function(score) {
  score <- as.matrix(score)

  row_max <- apply(score, 1L, max)

  exp_score <- exp(score - row_max)
  exp_score / rowSums(exp_score)
}

.afmvn_update_r <- function(M, B, J, m_q, Psi_q, nu_q) {
  
  score <- .afmvn_branch_scores(M=M, B = B, J = J, m_q = m_q, Psi_q = Psi_q, nu_q = nu_q)
  r_new <- .afmvn_row_softmax(score)

  list(r = r_new, score = score)
}