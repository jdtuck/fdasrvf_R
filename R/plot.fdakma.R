#' @export
plot.fdakma <- function(out){
  K <- length(out$fn)
  colors <- rainbow(K, s = 0.75)
  matplot(out$time, out$templates, type="l", col=colors, lty=1)
  title(main="Cluster Mean Functions")

  matplot(out$time, out$fn[[1]], type="l", col=colors()[349])
  for (k in 2:K){
    matplot(out$time, out$fn[[k]], type="l", col=colors()[347-4*k], add=T)
  }
  for (k in 1:K){
    lines(out$time, out$templates[,k], col=colors[k])
  }
  title(main="Clustered Functions")

  matplot(out$time, out$qn[[1]], type="l", col=colors()[349])
  for (k in 2:K){
    matplot(out$time, out$qn[[k]], type="l", col=colors()[347-4*k], add=T)
  }
  for (k in 1:K){
    lines(out$time, out$templates.q[,k], col=colors[k])
  }
  title(main="Clustered SRSF")
}