#' @export
plot.vfpca <- function(x, ...) {
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

  dims <- dim(x$q_pca)
  num.plot <- ceiling(dims[3] / 3)
  time <- x$time
  cnt <- 1
  stds <- x$stds
  idx <- which(stds == 0)
  for (ii in 1:num.plot) {
    graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE))

    graphics::matplot(
      time,
      x$q_pca[1:(dims[1] - 1), , cnt],
      type = "l",
      col = colors,
      xlab = "time",
      ylab = "aplitude"
    )
    graphics::lines(time, x$q_pca[1:(dims[1] - 1), idx, cnt], col = "black")
    graphics::title(main = sprintf("q domain: PD %d", cnt))

    if (dims[3] >= cnt + 1) {
      graphics::matplot(
        time,
        x$q_pca[1:(dims[1] - 1), , cnt + 1],
        type = "l",
        col = colors,
        xlab = "time",
        ylab = "aplitude"
      )
      graphics::lines(time, x$q_pca[1:(dims[1] - 1), idx, cnt + 1], col =
                        "black")
      graphics::title(main = sprintf("q domain: PD %d", cnt + 1))
    } else {
      graphics::plot.new()
    }
    if (dims[3] >= cnt + 2) {
      graphics::matplot(
        time,
        x$q_pca[1:(dims[1] - 1), , cnt + 2],
        type = "l",
        col = colors,
        xlab = "time",
        ylab = "aplitude"
      )
      graphics::lines(time, x$q_pca[1:(dims[1] - 1), idx, cnt + 2], col =
                        "black")
      graphics::title(main = sprintf("q domain: PD %d", cnt + 2))
    } else {
      graphics::plot.new()
    }

    graphics::matplot(
      time,
      x$f_pca[1:(dims[1] - 1), , cnt],
      type = "l",
      col = colors,
      xlab = "time",
      ylab = "aplitude"
    )
    graphics::lines(time, x$f_pca[1:(dims[1] - 1), idx, cnt], col = "black")
    graphics::title(main = sprintf("f domain: PD %d", cnt))
    if (dims[3] >= cnt + 1) {
      graphics::matplot(
        time,
        x$f_pca[1:(dims[1] - 1), , cnt + 1],
        type = "l",
        col = colors,
        xlab = "time",
        ylab = "aplitude"
      )
      graphics::lines(time, x$f_pca[1:(dims[1] - 1), idx, cnt + 1], col =
                        "black")
      graphics::title(main = sprintf("f domain: PD %d", cnt + 1), )
    } else {
      graphics::plot.new()
    }
    if (dims[3] >= cnt + 2) {
      graphics::matplot(
        time,
        x$f_pca[1:(dims[1] - 1), , cnt + 2],
        type = "l",
        col = colors,
        xlab = "time",
        ylab = "aplitude"
      )
      graphics::lines(time, x$f_pca[1:(dims[1] - 1), idx, cnt + 2], col =
                        "black")
      graphics::title(main = sprintf("f domain: PD %d", cnt + 2))
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
