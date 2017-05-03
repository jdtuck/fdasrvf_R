#' @export
plot.vfpca <- function(out){
  layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
  dims = dim(out$q_pca)
  time <- out$time
  matplot(time,out$q_pca[1:(dims[1]-1),,1],type="l")
  title(main="q domain: PD 1")
  matplot(time,out$q_pca[1:(dims[1]-1),,2],type="l")
  title(main="q domain: PD 2")
  matplot(time,out$q_pca[1:(dims[1]-1),,3],type="l")
  title(main="q domain: PD 3")
  matplot(time,out$f_pca[,,1],type="l")
  title(main="f domain: PD 1")
  matplot(time,out$f_pca[,,2],type="l")
  title(main="f domain: PD 2")
  matplot(time,out$f_pca[,,3],type="l")
  title(main="f domain: PD 3")
  layout(1)
  cumm_coef = 100*cumsum(out$latent)/sum(out$latent)
  plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")

}
