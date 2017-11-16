#' @export
plot.hfpca <- function(x, ...){
  dims <- dim(x$gam_pca)
  num.plot <- ceiling(dims[3]/3)
  TT <- dim(x$gam_pca)[2]
  time <- seq(0,1,len=TT)
  cnt <- 1
  for (ii in 1:num.plot){
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))

    matplot(time,t(x$gam_pca[,,cnt]),type="l")
    title(main=sprintf("gam: PD %d", cnt))
    if (dims[3] >= cnt + 1){
      matplot(time,t(x$gam_pca[,,cnt+1]),type="l")
      title(main=sprintf("gam: PD %d", cnt+1))
    } else {
      plot.new()
    }
    if (dims[3] >= cnt + 2){
      matplot(time,t(x$gam_pca[,,cnt+2]),type="l")
      title(main=sprintf("gam: PD %d", cnt+2))
    }else {
      plot.new()
    }

    cnt <- cnt + 3
  }

  layout(1)
  cumm_coef = 100*cumsum(x$latent)/sum(x$latent)
  plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")
}
