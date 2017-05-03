#' @export
plot.hfpca <- function(out){
  layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
  TT <- dim(out$gam_pca)[2]
  matplot(seq(0,1,len=TT),t(out$gam_pca[,,1]),type="l",xlab = "t",ylab = "t")
  title(main="PD 1")
  matplot(seq(0,1,len=TT),t(out$gam_pca[,,2]),type="l",xlab = "t",ylab = "t")
  title(main="PD 2")
  matplot(seq(0,1,len=TT),t(out$gam_pca[,,3]),type="l",xlab = "t",ylab = "t")
  title(main="PD 3")
  layout(1)
  cumm_coef = 100*cumsum(out$latent)/sum(out$latent)
  plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")
}
