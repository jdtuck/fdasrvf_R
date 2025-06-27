#' @export
plot.hfpns <- function(x, ...) {
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

  dims <- dim(x$gam_pns)
  num.plot <- ceiling(dims[3] / 3)
  TT <- dim(x$gam_pns)[2]
  time <- seq(0, 1, len = TT)
  cnt <- 1
  stds <- x$stds
  idx <- which(stds == 0)
  for (ii in 1:num.plot) {
    graphics::layout(matrix(c(1, 2, 3), 1, 3, byrow = TRUE))

    graphics::matplot(
      time,
      t(x$gam_pns[, , cnt]),
      type = "l",
      col = colors,
      xlab = "time",
      ylab = "warped time"
    )
    graphics::lines(time, t(x$gam_pns[idx, , cnt]), col = "black")
    graphics::title(main = sprintf("gam: PD %d", cnt))
    if (dims[3] >= cnt + 1) {
      graphics::matplot(
        time,
        t(x$gam_pns[, , cnt + 1]),
        type = "l",
        col = colors,
        xlab = "time",
        ylab = "warped time"
      )
      graphics::lines(time, t(x$gam_pns[idx, , cnt + 1]), col = "black")
      graphics::title(main = sprintf("gam: PD %d", cnt + 1))
    } else {
      graphics::plot.new()
    }
    if (dims[3] >= cnt + 2) {
      graphics::matplot(
        time,
        t(x$gam_pns[, , cnt + 2]),
        type = "l",
        col = colors,
        xlab = "time",
        ylab = "warped time"
      )
      graphics::lines(time, t(x$gam_pns[idx, , cnt + 2]), col = "black")
      graphics::title(main = sprintf("gam: PD %d", cnt + 2))
    } else {
      graphics::plot.new()
    }

    cnt <- cnt + 3
  }

  graphics::layout(1)
  cumm_coef = 100 * x$cumvar
  plot(
    cumm_coef,
    type = "l",
    col = colors[1],
    main = "Coefficient Cumulative Percentage",
    ylab = "Percentage"
  )
}
