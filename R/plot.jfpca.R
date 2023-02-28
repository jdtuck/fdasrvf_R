#' @export
plot.jfpca <- function(x, ...){
    dims <- dim(x$q_pca)
    num.plot <- ceiling(dims[3]/3)
    cnt <- 1
    for (ii in 1:num.plot){
      graphics::layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))

        time <- x$time
        graphics::matplot(time,x$q_pca[,,cnt],type="l")
        graphics::title(main=sprintf("q domain: PD %d", cnt))
        if (dims[3] >= cnt + 1){
          graphics::matplot(time,x$q_pca[,,cnt+1],type="l")
          graphics::title(main=sprintf("q domain: PD %d", cnt+1))
        } else {
          graphics::plot.new()
        }
        if (dims[3] >= cnt + 2){
          graphics::matplot(time,x$q_pca[,,cnt+2],type="l")
          graphics::title(main=sprintf("q domain: PD %d", cnt+2))
        }else {
          graphics::plot.new()
        }

        graphics::matplot(time,x$f_pca[,,cnt],type="l")
        graphics::title(main=sprintf("f domain: PD %d", cnt))
        if (dims[3] >= cnt + 1){
          graphics::matplot(time,x$f_pca[,,cnt+1],type="l")
          graphics::title(main=sprintf("f domain: PD %d", cnt+1))
        }else {
          graphics::plot.new()
        }
        if (dims[3] >= cnt + 2){
          graphics::matplot(time,x$f_pca[,,cnt+2],type="l")
          graphics::title(main=sprintf("f domain: PD %d", cnt+2))
        }

        cnt <- cnt + 3
    }

    graphics::layout(1)
    cumm_coef <- 100*cumsum(x$latent)/sum(x$latent)
    plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")

}
