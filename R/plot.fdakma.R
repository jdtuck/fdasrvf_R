#' @export
plot.fdakma <- function(x, ...){
  K <- dim(x$templates)[2]
  num.plot <- ceiling(K/6)
  colors <- grDevices::rainbow(K, s = 0.75)
  graphics::matplot(x$time, x$templates, type="l", col=colors, lty=1)
  graphics::title(main="Cluster Mean Functions")

  for (k in 1:num.plot){
    graphics::layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
    for (n in ((k-1)*6+1):min(K,(k)*6)){
      graphics::matplot(x$time, x$fn[[n]], type="l", col=colors()[349], lty=1,
              xlab="Time",ylab="")
      graphics::lines(x$time, x$templates[,n], col=colors[n])
      graphics::title(main=sprintf("Cluster f: %d",n))
    }
  }

  for (k in 1:num.plot){
    graphics::layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
    for (n in ((k-1)*6+1):min(K,(k)*6)){
      graphics::matplot(x$time, x$qn[[n]], type="l", col=colors()[349], lty=1,
              xlab="Time",ylab="")
      graphics::lines(x$time, x$templates.q[,n], col=colors[n])
      graphics::title(main=sprintf("Cluster q: %d",n))
    }
  }

}
