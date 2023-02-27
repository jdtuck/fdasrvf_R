#' @export
plot.ampbox <- function(x, ...){
  fmedian <- x$fmedian
  maxx <- x$maxx
  minn <- x$minn
  Q1 <- x$Q1
  Q1a <- x$Q1a
  Q3 <- x$Q3
  Q3a <- x$Q3a
  time <- x$time
  M <- length(fmedian)
  ymin <- min(c(min(fmedian),min(Q1),min(Q3),min(maxx),min(minn)))
  ymax <- max(c(max(fmedian),max(Q1),max(Q3),max(maxx),max(minn)))
  plot(time, fmedian, col="black",xlab="Time",main="Amplitude Boxplot", type="l", ylim=c(ymin, ymax))
  graphics::lines(time, Q1, col="blue")
  graphics::lines(time, Q3, col="blue")
  graphics::lines(time, Q1a, col="green")
  graphics::lines(time, Q3a, col="green")
  graphics::lines(time, maxx, col="red")
  graphics::lines(time, minn, col="red")

  s <- seq(0,1,length.out=100)
  Fs2 <- matrix(0,length(time), 595)
  Fs2[,1] <- (1-s[1]) * minn + s[1] * Q1
  for (j in 2:100){
    Fs2[,j] <- (1-s[j]) * minn + s[j] * Q1a
    Fs2[,99+j] <- (1-s[j]) * Q1a + s[j] * Q1
    Fs2[,198+j] <- (1-s[j]) * Q1 + s[j] * fmedian
    Fs2[,297+j] <- (1-s[j]) * fmedian + s[j] * Q3
    Fs2[,396+j] <- (1-s[j]) * Q3 + s[j] * Q3a
    Fs2[,495+j] <- (1-s[j]) * Q3a + s[j] * maxx
  }
  d1<-sqrt(trapz(time,(x$qmedian-x$Q1_q)^2))
  d1a<-sqrt(trapz(time,(x$Q1_q-x$Q1a_q)^2))
  dl<-sqrt(trapz(time,(x$Q1a_q-x$min_q)^2))
  d3<-sqrt(trapz(time,(x$qmedian-x$Q3_q)^2))
  d3a<-sqrt(trapz(time,(x$Q3_q-x$Q3a_q)^2))
  du<-sqrt(trapz(time,(x$Q3a_q-x$max_q)^2))
  part1<-seq(-d1-d1a-dl,-d1-d1a,length.out=100)
  part2<-seq(-d1-d1a,-d1,length.out=100)
  part3<-seq(-d1,0,length.out=100)
  part4<-seq(0,d3,length.out=100)
  part5<-seq(d3,d3+d3a,length.out=100)
  part6<-seq(d3+d3a,d3+d3a+du,length.out=100)
  allparts<-c(part1,part2[2:100],part3[2:100],part4[2:100],part5[2:100],part6[2:100])

  if (requireNamespace("plot3Drgl", quietly = TRUE)) {
    p <- plot3D::persp3D(x=time,y=allparts,z=Fs2,col= viridisLite::viridis(128),plot=F,main="Amplitude Surface Plot",ticktype="detailed",box=F)+
      plot3D::lines3D(x=time,y=rep(0,M),z=fmedian,col="black",lwd=6,add=T,plot=F)+
      plot3D::lines3D(x=time,y=rep(-d1,M),z=Q1,col="blue",lwd=6,add=T,plot=F)+
      plot3D::lines3D(x=time,y=rep(-d1-d1a,M),z=Q1a,col="green",lwd=6,add=T,plot=F)+
      plot3D::lines3D(x=time,y=rep(-d1-d1a-dl, M),z=minn,col="red",lwd=6,add=T,plot=F)+
      plot3D::lines3D(x=time,y=rep(d3, M),z=Q3,col="blue",lwd=6,add=T,plot=F)+
      plot3D::lines3D(x=time,y=rep(d3+d3a, M),z=Q3a,col="green",lwd=6,add=T,plot=F)+
      plot3D::lines3D(x=time,y=rep(d3+d3a+du, M),z=maxx,col="red",lwd=6,add=T,plot=F)
    plot3Drgl::plotrgl()
    rgl::par3d("windowRect"= c(0,0,640,640))
    rgl::grid3d(c("x", "y+", "z"))
    rgl::axes3d(c('x--',"y--",'z'))
    rgl::title3d(xlab="Time",ylab="Distance")
  } else {
    graphics::image(time, allparts, Fs2, main="Surface Plot", ylab="", col=viridisLite::viridis(128))
    graphics::lines(time, rep(0, M), col="black", lwd=1)
    graphics::lines(time, rep(-d1, M), col="blue", lwd=1)
    graphics::lines(time, rep(-d1-d1a, M), col="green", lwd=1)
    graphics::lines(time, rep(-d1-d1a-dl, M), col="red", lwd=1)
    graphics::lines(time, rep(d3, M), col="blue", lwd=1)
    graphics::lines(time, rep(d3+d3a, M), col="green", lwd=1)
    graphics::lines(time, rep(d3+d3a+du, M), col="red", lwd=1)
  }
}
