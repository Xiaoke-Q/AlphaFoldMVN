#' Visualize the data points of length 3 and boundaries.
#'
#' @param W A n-by-3 matrix.
#' @param alpha A numeric parameter of alpha transformation.
#' @param group A vector of length n, labeling the groups.
#' @param shift A vector of length n, labeling the shifts.
#' @param color A vector of length n, labeling the colors.
#' @param opacity A vector of length n, labeling the transparency.
#'
#' @returns A plotly object.
#'
#' @export
#' @examples visualize(W= W, alpha = 0.7, group = rep(1:3, each = nrow(W)))
visualize <- function(W = NULL, alpha, group = NULL, shift = NULL, color = NULL, opacity = NULL){
  if(!is.null(W)){
    span <- max(c(abs(W)), abs(1/alpha))
  }
  span <- abs(1/alpha)
  span <- span * 3
  xs <- seq(-span, span, length.out = 60)
  ys <- seq(-span, span, length.out = 60)
  zs  <- outer(xs, ys, function(x, y) -x - y)
  
  fig <- plotly::plot_ly(x = ~xs, y = ~ys, z = ~zs) |>
    plotly::add_surface(opacity = 0.6, showscale = FALSE,
                colorscale = list(c(0,1), c("lightblue","lightblue")),
                name = "x+y+z=0")
  
  y_line <- seq(-(1/alpha), (2/alpha), length.out = 100)
  z_line <- 1/alpha - y_line
  x_line <- rep(-1/alpha, length(y_line))
  
  fig <- fig |>
    plotly::add_trace(x = ~x_line, y = ~y_line, z = ~z_line,
              type = "scatter3d", mode = "lines",
              line = list(color = "#4F2C1D", width = 6))
  
  x_line2 <- seq(-(1/alpha), (2/alpha), length.out = 100)
  z_line2 <- 1/alpha - x_line2
  y_line2 <- rep(-1/alpha, length(x_line2))
  
  fig <- fig |>
    plotly::add_trace(x = ~x_line2, y = ~y_line2, z = ~z_line2,
              type = "scatter3d", mode = "lines",
              line = list(color = "#4F2C1D", width = 6))
  y_line3 <- seq(-(1/alpha), (2/alpha), length.out = 100)
  x_line3 <- 1/alpha - y_line3
  z_line3 <- rep(-1/alpha, length(y_line3))
  
  fig <- fig |>
    plotly::add_trace(x = ~x_line3, y = ~y_line3, z = ~z_line3,
              type = "scatter3d", mode = "lines",
              line = list(color = "#4F2C1D", width = 6))
  
  pal <- c("#F2A900","#385E9D","#7FB5E3", "#F28E2B", "#8C564B", "#E15759")

  if(!is.null(W)){
    if(is.null(group)){
      group = rep(1, nrow(W))
    }

    ng <- length(unique(group))
    ds <- if(ng == 1) c(0) else seq(-span/4, span/4, length.out = ng)
    os <- if(ng == 1) c(1) else seq(0, 1, length.out = ng)
    if(is.null(color)) color <- pal[group]
    if(is.null(shift)) shift <- ds[group]
    if(is.null(opacity)) opacity <- os[group]
    

    W <- W + shift
    fig <- fig |>
      plotly::add_markers(
        x = W[,1], 
        y = W[,2], 
        z = W[,3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 2, 
          color = color,
          opacity = opacity
        ))
  }

  fig

}

#' Visualize observed data X
#'
#' @param X A matrix
#' @param group Group of data
#' @param color Color
#'
#' @returns A plotly object
#'
#' @export
#' @examples visualize_x(X = X)
visualize_x <- function(X, group = NULL, color = NULL){

  xs <- seq(0, 1, length.out = 60)
  ys <- seq(0, 1, length.out = 60)
  zs  <- outer(xs, ys, function(x, y) 1 - x - y)
  zs[zs < 0] <- NA
  fig <- plotly::plot_ly(x = ~xs, y = ~ys, z = ~zs) |>
    plotly::add_surface(opacity = 0.6, showscale = FALSE,
                colorscale = list(c(0,1), c("lightblue","lightblue")),
                name = "x+y+z=0")
  
  y_line <- seq(0, 1, length.out = 100)
  z_line <- 1 - y_line
  x_line <- rep(0, length(y_line))
  
  fig <- fig |>
    plotly::add_trace(x = ~x_line, y = ~y_line, z = ~z_line,
              type = "scatter3d", mode = "lines",
              line = list(color = "#4F2C1D", width = 6))
  
  x_line2 <- seq(0, 1, length.out = 100)
  z_line2 <- 1 - x_line2
  y_line2 <- rep(0, length(x_line2))
  
  fig <- fig |>
    plotly::add_trace(x = ~x_line2, y = ~y_line2, z = ~z_line2,
              type = "scatter3d", mode = "lines",
              line = list(color = "#4F2C1D", width = 6))
  y_line3 <- seq(0, 1, length.out = 100)
  x_line3 <- 1 - y_line3
  z_line3 <- rep(0, length(y_line3))
  
  fig <- fig |>
    plotly::add_trace(x = ~x_line3, y = ~y_line3, z = ~z_line3,
              type = "scatter3d", mode = "lines",
              line = list(color = "#4F2C1D", width = 6))
  
  pal <- c("#F2A900","#385E9D","#7FB5E3", "#F28E2B", "#8C564B", "#E15759")
    if(is.null(group)){
      group = rep(1, nrow(X))
    }

    ng <- length(unique(group))
    if(is.null(color)) color <- pal[group]
    

  
  fig <- fig |>
    plotly::add_markers(x = X[,1], y = X[,2], z = X[,3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 2, 
          color = color
        ))
  
  fig
}
