#' Evaluate the variational eta objective on a grid
#'
#' Computes the eta-block variational objective over a grid of `m_eta` and
#' `log_s_eta` values. This is mainly intended for diagnostic visualization
#' before running the eta update optimizer.
#'
#' @param X An `N` by `D` compositional data matrix.
#' @param r An `N` by 2 matrix of branch responsibilities.
#' @param m_q Current variational mean vector for `mu`.
#' @param Psi_q Current Wishart scale matrix for `Lambda`.
#' @param nu_q Current Wishart degrees of freedom.
#' @param m_grid Numeric grid for `m_eta`.
#' @param log_s_grid Numeric grid for `log(s_eta)`.
#' @param a0 Prior mean for `eta`.
#' @param s0_sq Prior variance for `eta`.
#' @param n_quad Number of Gauss-Hermite quadrature nodes.
#' @param eps_alpha Threshold for treating `alpha` as zero.
#'
#' @return An object of class `"afmvn_eta_surface"`.
#' @export
eta_objective_surface <- function(X, r, m_q, Psi_q, nu_q,
                                  m_grid = seq(-4, 4, length.out = 80L),
                                  log_s_grid = seq(log(0.03), log(5), length.out = 80L),
                                  a0 = 0, s0_sq = 10,
                                  n_quad = 20L, eps_alpha = 1e-8) {
  X <- as.matrix(X)
  r <- as.matrix(r)
  Psi_q <- as.matrix(Psi_q)

  L <- matrix(NA_real_, nrow = length(log_s_grid), ncol = length(m_grid))

  for (j in seq_along(m_grid)) {
    for (k in seq_along(log_s_grid)) {
      L[k, j] <- .afmvn_eta_objective(
        X = X,
        m_eta = m_grid[j],
        s2_eta = exp(2 * log_s_grid[k]),
        r = r,
        m_q = m_q,
        Psi_q = Psi_q,
        nu_q = nu_q,
        a0 = a0,
        s0_sq = s0_sq,
        n_quad = n_quad,
        eps_alpha = eps_alpha
      )
    }
  }

  idx <- which(L == max(L, na.rm = TRUE), arr.ind = TRUE)

  best <- list(
    m_eta = m_grid[idx[1L, "col"]],
    log_s_eta = log_s_grid[idx[1L, "row"]],
    s_eta = exp(log_s_grid[idx[1L, "row"]]),
    s2_eta = exp(2 * log_s_grid[idx[1L, "row"]]),
    value = L[idx[1L, "row"], idx[1L, "col"]]
  )

  out <- list(
    L = L,
    L_relative = L - max(L, na.rm = TRUE),
    m_grid = m_grid,
    log_s_grid = log_s_grid,
    s_grid = exp(log_s_grid),
    s2_grid = exp(2 * log_s_grid),
    best = best,
    call = match.call()
  )

  class(out) <- "afmvn_eta_surface"
  out
}

#' Plot an eta objective surface
#'
#' @param x An object returned by [eta_objective_surface()].
#' @param relative Logical. If `TRUE`, plot the objective relative to its grid maximum.
#' @param nlevels Number of contour levels.
#' @param mark_best Logical. If `TRUE`, mark the best grid point.
#' @param image Logical. If `TRUE`, draw an image heatmap before contours.
#' @param ... Additional argument.
#'
#' @return Invisibly returns `x`.
#' @method plot afmvn_eta_surface
#' @export
plot.afmvn_eta_surface <- function(x,
                                   relative = TRUE,
                                   nlevels = 20L,
                                   mark_best = TRUE,
                                   image = TRUE,
                                   ...) {
  if (!inherits(x, "afmvn_eta_surface")) {
    stop("x must be an 'afmvn_eta_surface' object.", call. = FALSE)
  }

  Z <- if (isTRUE(relative)) x$L_relative else x$L

  if (isTRUE(image)) {
    graphics::image(
      x = x$m_grid,
      y = x$log_s_grid,
      z = t(Z),
      xlab = "m_eta",
      ylab = "log(s_eta)",
      main = "Eta objective"
    )

    graphics::contour(
      x = x$m_grid,
      y = x$log_s_grid,
      z = t(Z),
      add = TRUE,
      nlevels = nlevels,
      ...
    )
  }else{
    graphics::contour(
      x = x$m_grid,
      y = x$log_s_grid,
      z = t(Z),
      xlab = "m_eta",
      ylab = "log(s_eta)",
       main = "Eta objective",
      nlevels = nlevels,
      ...
    )
  }

  if (isTRUE(mark_best)) {
    graphics::points(
      x$best$m_eta,
      x$best$log_s_eta,
      pch = 19
    )
  }

  invisible(x)
}