#' @export
plot.curve_pca <- function(x, ...) {
  colors = c(
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#B3B3B3"
  )

  dims <- dim(x$pd)
  num.plot <- ceiling(dims[1] / 3)
  stds <- x$stds
  idx <- which(stds == 0)
  cnt = 1
  for (ii in 1:num.plot) {
    graphics::layout(matrix(c(1, 2, 3), 1, 3, byrow = TRUE))

    plot_curve(x$pd[cnt, idx][[1]], col="black")
    for (jj in 1:dims[2]){
      if (jj == idx){
        next
      }
      plot_curve(x$pd[cnt, jj][[1]], add=TRUE, col=colors[jj])
    }
    graphics::title(main = sprintf("PD %d", cnt))
    if (dims[1] >= cnt + 1) {
      plot_curve(x$pd[cnt + 1, idx][[1]], col="black")
      for (jj in 1:dims[2]){
        if (jj == idx){
          next
        }
        plot_curve(x$pd[cnt+1, jj][[1]], add=TRUE, col=colors[jj])
      }
      graphics::title(main = sprintf("PD %d", cnt + 1))
    } else {
      graphics::plot.new()
    }
    if (dims[1] >= cnt + 2) {
      plot_curve(x$pd[cnt + 1, idx][[1]], col="black")
      for (jj in 1:dims[2]){
        if (jj == idx){
          next
        }
        plot_curve(x$pd[cnt + 2, jj][[1]], add=TRUE, col=colors[jj])
      }
      graphics::title(main = sprintf("PD %d", cnt + 2))
    } else {
      graphics::plot.new()
    }

    cnt <- cnt + 3
  }

  graphics::layout(1)
  cumm_coef = 100 * cumsum(x$latent) / sum(x$latent)
  plot(
    cumm_coef,
    type = "l",
    col = colors[1],
    main = "Coefficient Cumulative Percentage",
    ylab = "Percentage"
  )
}
