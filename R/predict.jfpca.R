#' @export
predict.jfpca <- function(x, f, ...){
  q1 = f_to_srvf(f, x$warp_data$time)
  M = length(x$warp_data$time)
  N = ncol(f)
  gam = matrix(0,M,N)
  fn = matrix(0,M,N)
  qn = matrix(0,M,N)
  for (ii in 1:N){
    gam[,ii] = optimum.reparam(x$warp_data$mqn, x$warp_data$time, q1[,ii],
                               x$warp_data$time, method=x$warp_data$call$optim_method)
    fn[,ii] = warp_f_gamma(f[,ii], x$warp_data$time, gam[,ii])
    qn[,ii] = f_to_srvf(fn[,ii], x$warp_data$time)
  }

  m_new <- sign(fn[x$id,])*sqrt(abs(fn[x$id,]))  # scaled version
  qn1 <- rbind(qn,m_new)

  no = ncol(x$U)
  psi = matrix(0,M,N)
  vec = matrix(0,M,N)
  binsize <- mean(diff(x$warp_data$time))
  for (i in 1:N){
    psi[,i] = sqrt(gradient(gam[,i],binsize))
    vec[,i] <- inv_exp_map(x$mu_psi, psi[,i])
  }

  g <- rbind(qn1,x$C*vec)
  a <- matrix(0,N,no)
  for (i in 1:N){
    for (j in 1:no){
      a[i,j] <- (g[,i]-x$mu_g)%*%x$U[,j]
    }
  }

  a
}
