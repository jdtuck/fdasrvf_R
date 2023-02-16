#' Vertical Functional Principal Component Analysis
#'
#' This function calculates vertical functional principal component analysis
#' on aligned data
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param no number of principal components to extract
#' @param id point to use for f(0) (default = midpoint)
#' @param ci geodesic standard deviations (default = c(-1,0,1))
#' @param showplot show plots of principal directions (default = T)
#' @return Returns a vfpca object containing \item{q_pca}{srvf principal directions}
#' \item{f_pca}{f principal directions}
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' \item{id}{point used for f(0)}
#' @keywords srvf alignment
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_warp")
#' vfpca = vertFPCA(simu_warp,no = 3)
vertFPCA <- function(warp_data,no,id=round(length(warp_data$time)/2),ci=c(-1,0,1),showplot = TRUE){
    # Parameters
    fn <- warp_data$fn
    time <- warp_data$time
    qn <- warp_data$qn
    NP = 1:no  # number of principal components
    Nstd = length(ci)

    # FPCA
    mq_new = rowMeans(qn)
    m_new = sign(fn[id,])*sqrt(abs(fn[id,]))  # scaled version
    mqn = c(mq_new,mean(m_new))
    K = cov(t(rbind(qn,m_new))) #out$sigma

    out = svd(K)
    s = out$d
    stdS = sqrt(s)
    U = out$u

    # compute the PCA in the q domain
    q_pca = array(0,dim=c((length(mq_new)+1),Nstd,no))
    for (k in NP){
        for (i in 1:Nstd){
            q_pca[,i,k] = mqn + ci[i]*stdS[k]*U[,k]
        }
    }

    # compute the correspondence to the original function domain
    f_pca = array(0,dim=c((length(mq_new)),Nstd,no))
    for (k in NP){
        for (i in 1:Nstd){
            if (id == 1){
              f_pca[,i,k] <- cumtrapz(time,q_pca[1:(dim(q_pca)[1]-1),i,k]*
                                      abs(q_pca[1:(dim(q_pca)[1]-1),i,k]))+(sign(q_pca[dim(q_pca)[1],i,k])*(q_pca[dim(q_pca)[1],i,k]^2))

            } else {
              f_pca[,i,k] <- cumtrapzmid(time,q_pca[1:(dim(q_pca)[1]-1),i,k]*
                                          abs(q_pca[1:(dim(q_pca)[1]-1),i,k]),sign(q_pca[dim(q_pca)[1],i,k])*
                                          (q_pca[dim(q_pca)[1],i,k]^2), id)
            }
        }
        fbar = rowMeans(fn)
        fsbar = rowMeans(f_pca[,,k])
        err = kronecker(matrix(1,1,Nstd),fbar-fsbar)
        f_pca[,,k] = f_pca[,,k] + err
    }

    N2 = dim(qn)[2]
    c = matrix(0,N2,no)
    for (k in NP){
        for (i in 1:N2){
            c[i,k] = sum((c(qn[,i],m_new[i])-mqn)*U[,k])
        }
    }

    vfpca <- list()
    vfpca$q_pca <- q_pca
    vfpca$f_pca <- f_pca
    vfpca$latent <- s[NP]
    vfpca$coef <- c[,NP]
    vfpca$U <- U[,NP]
    vfpca$id <- id
    vfpca$mqn <- mqn
    vfpca$time <- time

    class(vfpca) <- "vfpca"

    if (showplot){
        plot(vfpca)
    }

    return(vfpca)
}
