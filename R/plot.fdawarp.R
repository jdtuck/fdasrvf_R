#' @export
plot.fdawarp <- function(x, ...){
  M <- nrow(x$fn)

  mean_f0 <- rowMeans(x$f0);
  std_f0 <- apply(x$f0, 1, stats::sd)
  mean_fn <- rowMeans(x$fn)
  std_fn <- apply(x$fn, 1, stats::sd)
  graphics::matplot((0:(M-1))/(M-1),x$gam,type="l",main="Warping functions",xlab="Time")

  graphics::matplot(x$time,x$fn,type="l",main=bquote(paste("Warped Data ",lambda ==
      .(x$lambda))))

  graphics::matplot(x$time,cbind(mean_f0,mean_f0+std_f0,mean_f0-std_f0),type="l",lty=1,
                  col=c("blue","red","green"),
                  ylab="",main=bquote(paste("Original Data: ", Mean %+-% STD)))
  graphics::legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
               col=c('blue','red','green'),lty=1)

  graphics::matplot(x$time,cbind(mean_fn,mean_fn+std_fn,mean_fn-std_fn),type="l",lty=1,
                  col=c("blue","red","green"),
                  ylab="",main=bquote(paste("Warped Data: ",lambda == .(x$lambda),": ",
                                                                      Mean %+-% STD)))
  graphics::legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
               col=c('blue','red','green'),lty=1)


  if (x$method=="mean"){
    plot(x$time,x$fmean,type="l",col="green",main=bquote(paste(f[mean]," ", lambda == .(x$lambda))))
  } else {
    plot(x$time,x$fmean,type="l",col="green",main=bquote(paste(f[median]," ", lambda == .(x$lambda))))
  }

}
