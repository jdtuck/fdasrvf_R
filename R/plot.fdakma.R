#' @export
plot.fdakma <- function(x, ...){
  K <- length(x$fn)
  colors <- rainbow(K, s = 0.75)
  matplot(x$time, x$templates, type="l", col=colors, lty=1)
  title(main="Cluster Mean Functions")

  matplot(x$time, x$fn[[1]], type="l", col=colors()[349])
  for (k in 2:K){
    matplot(x$time, x$fn[[k]], type="l", col=colors()[347-4*k], add=T)
  }
  for (k in 1:K){
    lines(x$time, x$templates[,k], col=colors[k])
  }
  title(main="Clustered Functions")

  matplot(x$time, x$qn[[1]], type="l", col=colors()[349])
  for (k in 2:K){
    matplot(x$time, x$qn[[k]], type="l", col=colors()[347-4*k], add=T)
  }
  for (k in 1:K){
    lines(x$time, x$templates.q[,k], col=colors[k])
  }
  title(main="Clustered SRSF")
}