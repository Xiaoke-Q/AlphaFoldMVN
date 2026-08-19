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

.afmvn_quad_from_mbj <- function(mbj, quad_ctx) {
  M <- mbj$M
  B <- mbj$B

  N <- dim(M)[1L]
  G <- dim(M)[2L]
  p <- dim(M)[3L]

  M_flat <- matrix(M, nrow = N * G, ncol = p)
  B_flat <- matrix(B, nrow = N * G, ncol = p * p)

  quad <- drop(B_flat %*% quad_ctx$vec_Psi_t) -
    2 * drop(M_flat %*% quad_ctx$Psi_m) +
    quad_ctx$c_q

  matrix(quad, nrow = N, ncol = G)
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

.afmvn_update_r <- function(mbj, nu_q, quad = NULL,
                            quad_ctx = NULL, m_q = NULL, Psi_q = NULL) {
  if (is.null(quad)) {
    if (is.null(quad_ctx)) {
      if (is.null(m_q) || is.null(Psi_q)) {
        stop("Provide either quad, quad_ctx, or both m_q and Psi_q.", call. = FALSE)
      }
      quad_ctx <- .afmvn_quad_context(m_q = m_q, Psi_q = Psi_q)
    }

    quad <- .afmvn_quad_from_mbj(mbj = mbj, quad_ctx = quad_ctx)
  }

  score <- .afmvn_branch_scores(mbj = mbj, quad = quad, nu_q = nu_q)
  r_new <- .afmvn_row_softmax(score)

  list(r = r_new, score = score, quad = quad)
}