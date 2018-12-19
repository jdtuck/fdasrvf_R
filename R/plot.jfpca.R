#' @export
plot.jfpca <- function(out){
    dims <- dim(out$q_pca)
    num.plot <- ceiling(dims[3]/3)
    cnt <- 1
    for (ii in 1:num.plot){
        layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))

        time <- out$time
        matplot(time,out$q_pca[,,cnt],type="l")
        title(main=sprintf("q domain: PD %d", cnt))
        if (dims[3] >= cnt + 1){
            matplot(time,out$q_pca[,,cnt+1],type="l")
            title(main=sprintf("q domain: PD %d", cnt+1))
        } else {
            plot.new()
        }
        if (dims[3] >= cnt + 2){
            matplot(time,out$q_pca[,,cnt+2],type="l")
            title(main=sprintf("q domain: PD %d", cnt+2))
        }else {
            plot.new()
        }

        matplot(time,out$f_pca[,,cnt],type="l")
        title(main=sprintf("f domain: PD %d", cnt))
        if (dims[3] >= cnt + 1){
            matplot(time,out$f_pca[,,cnt+1],type="l")
            title(main=sprintf("f domain: PD %d", cnt+1))
        }else {
            plot.new()
        }
        if (dims[3] >= cnt + 2){
            matplot(time,out$f_pca[,,cnt+2],type="l")
            title(main=sprintf("f domain: PD %d", cnt+2))
        }

        cnt <- cnt + 3
    }

    layout(1)
    cumm_coef <- 100*cumsum(out$latent)/sum(out$latent)
    plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")

}
