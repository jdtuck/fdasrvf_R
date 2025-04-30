#' @importFrom graphics plot
#' @export
plot.fdacurve <- function(x, ...){

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

  N <- dim(x$beta)[3]

  plot_curve(x$beta[,,1], col=colors[1], main="Original Curves")
  for (i in 2:N){
    k = i %% length(colors)
    plot_curve(x$beta[,,i], add=TRUE, col=colors[k])
  }

  idx = which.max(x$betan[1,1,])
  plot_curve(x$betan[,,idx], col=colors[1], main="Aligned Curves")
  for (i in 1:N){
    k = i %% length(colors)
    plot_curve(x$betan[,,i], add=TRUE, col=colors[k])
  }

  M = dim(x$beta)[2]
  graphics::matplot(
    x = (0:(M - 1)) / (M - 1),
    y = x$gam,
    type = "l",
    col = colors,
    main = "Warping functions",
    xlab = "Time",
    ylab = "Warped Time"
  )

}
