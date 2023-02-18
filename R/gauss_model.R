#' Gaussian model of functional data
#'
#' This function models the functional data using a Gaussian model extracted from
#' the principal components of the srvfs
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param n number of random samples (n = 1)
#' @param sort_samples sort samples (default = F)
#' @return Returns a fdawarp object containing \item{fs}{random aligned samples}
#' \item{gams}{random warping function samples}
#' \item{ft}{random function samples}
#' @keywords pca
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' out1 <- gauss_model(simu_warp, n = 10)
gauss_model <- function(warp_data,n = 1,sort_samples = FALSE){
    fn <- warp_data$fn
    time <- warp_data$time
    qn <- warp_data$qn
    gam <- warp_data$gam
    # Parameters
    no = 3
    eps = .Machine$double.eps
    binsize = mean(diff(time))
    M = length(time)

    # compute mean and covariance in q-domain
    mq_new = rowMeans(qn)
    id = round(length(time)/2)
    m_new = sign(fn[id,])*sqrt(abs(fn[id,]))  # scaled version
    mqn = c(mq_new,mean(m_new))
    C = cov(t(rbind(qn,m_new)))

    q_s = rmvnorm(n,mean=mqn,sigma=C,method="svd")
    q_s = t(q_s)
    end = dim(q_s)[1]

    # compute the correspondence to the original function domain
    fs = matrix(0,M,n)
    for (k in 1:n){
        fs[,k] = cumtrapzmid(time,q_s[1:(end-1),k]*abs(q_s[1:(end-1),k]),sign(q_s[end,k])*(q_s[end,k]^2),id)
    }
    fbar = rowMeans(fn)
    fsbar = rowMeans(fs)
    err = kronecker(matrix(1,1,n),fbar-fsbar)
    fs = fs + err

    # random warping generation
    rgam = randomGamma(gam,n)
    gams = matrix(0,n,M)
    for (k in 1:n){
        gams[k,] = invertGamma(rgam[k,])
    }
    gams = t(gams)

    # sort functions and warpings
    if (sort_samples == T){
        mx = apply(fs,2, max)
        out_sort = sort(mx,index.return=TRUE)
        seq1 = out_sort$ix

        # compute the psi-function
        fy = gradient(t(rgam),binsize)
        psi = fy/sqrt(abs(fy)+eps)
        psi = t(psi)
        ip = rep(0,n)
        len = rep(0,n)
        for (i in 1:n){
            ip[i] = rep(1,M)%*%psi[i,]/M;
            len[i] = acos(rep(1,M)%*%psi[i,]/M)
        }
        out_sort = sort(len,index.return=TRUE)
        seq2 = out_sort$ix

        # combine x-variability and y-variability
        ft = matrix(0,M,n)
        for (k in 1:n){
            tmp = approx((0:(M-1))/(M-1),fs[,seq1[k]],xout = gams[,seq2[k]])
            ft[,k] = tmp$y
            while (is.na(ft[,k])){
                rgam2 = randomGamma(gam,1)
                tmp = approx((0:(M-1))/(M-1),fs[,seq1[k]],xout = invertGamma(rgam2))
                ft[,k] = tmp$y
            }
        }
    }else
    {
        # combine x-variability and y-variability
        ft = matrix(0,M,n)
        for (k in 1:n){
            tmp = approx((0:(M-1))/(M-1),fs[,k],xout = gams[,k])
            ft[,k] = tmp$y
            while (is.na(ft[,k])[1]){
                rgam2 = randomGamma(gam,1)
                tmp = approx((0:(M-1))/(M-1),fs[,k],xout = invertGamma(rgam2))
                ft[,k] = tmp$y
            }
        }
    }

    warp_data$fs = fs
    warp_data$gams = rgam
    warp_data$ft = ft
    warp_data$qs = q_s[1:(end-1),]
    warp_data$rsamps=T

    return(warp_data)
}
