#' @export
plot.vfpca <- function(x, ...){
  layx(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
  dims = dim(x$q_pca)
  time <- x$time
  matplot(time,x$q_pca[1:(dims[1]-1),,1],type="l")
  title(main="q domain: PD 1")
  matplot(time,x$q_pca[1:(dims[1]-1),,2],type="l")
  title(main="q domain: PD 2")
  matplot(time,x$q_pca[1:(dims[1]-1),,3],type="l")
  title(main="q domain: PD 3")
  matplot(time,x$f_pca[,,1],type="l")
  title(main="f domain: PD 1")
  matplot(time,x$f_pca[,,2],type="l")
  title(main="f domain: PD 2")
  matplot(time,x$f_pca[,,3],type="l")
  title(main="f domain: PD 3")
  layx(1)
  cumm_coef = 100*cumsum(x$latent)/sum(x$latent)
  plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")

}
