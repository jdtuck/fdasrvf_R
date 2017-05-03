#' @export
plot.fdawarp <- function(out){
  M <- nrow(out$fn)

  mean_f0 <- rowMeans(out$f0);
  std_f0 <- apply(out$f0, 1, sd)
  mean_fn <- rowMeans(out$fn)
  std_fn <- apply(out$fn, 1, sd)
  matplot((0:(M-1))/(M-1),out$gam,type="l",main="Warping functions",xlab="Time")

  matplot(out$time,out$fn,type="l",main=bquote(paste("Warped Data ",lambda ==
      .(out$lambda))))

  matplot(out$time,cbind(mean_f0,mean_f0+std_f0,mean_f0-std_f0),type="l",lty=1,
                  col=c("blue","red","green"),
                  ylab="",main=bquote(paste("Original Data: ", Mean %+-% STD)))
  legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
               col=c('blue','red','green'),lty=1)

  matplot(out$time,cbind(mean_fn,mean_fn+std_fn,mean_fn-std_fn),type="l",lty=1,
                  col=c("blue","red","green"),
                  ylab="",main=bquote(paste("Warped Data: ",lambda == .(out$lambda),": ",
                                                                      Mean %+-% STD)))
  legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
               col=c('blue','red','green'),lty=1)


  if (out$method=="mean"){
    plot(out$time,out$fmean,type="l",col="green",main=bquote(paste(f[mean]," ", lambda == .(out$lambda))))
  } else {
    plot(out$time,out$fmean,type="l",col="green",main=bquote(paste(f[median]," ", lambda == .(out$lambda))))
  }

}
