#' Joint Vertical and Horizontal Functional Principal Component Analysis
#'
#' This function calculates amplitude and phase joint  functional principal component
#' analysis on aligned data
#'
#' @param warp_data fdawarp objecet from \link{time_warping} of aligned data
#' @param no number of prinicpal components to extract
#' @param id integration point for f0 (default = midpoint)
#' @param C balance value (default = NULL)
#' @param showplot show plots of prinipal directions (default = T)
#' @return Returns a list containing \item{q_pca}{srvf principal directions}
#' \item{f_pca}{f principal directions}
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' \item{mu_psi}{mean psi function}
#' \item{mu_g}{mean g function}
#' \item{id}{point use for f(0)}
#' \item{C}{optimized phase amplitude ratio}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Jung, S. L. a. S. (2016). "Combined Analysis of Amplitude and Phase Variations in Functional Data."
#'        	arXiv:1603.01775 [stat.ME].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_warp")
#' data("simu_data")
#' jfpca = jointFPCA(simu_warp$fn,simu_data$time,simu_warp$qn,
#'                   simu_warp$gam, no = 3)
jointFPCA <- function(warp_data, no, id=round(length(warp_data$time)/2), C=NULL, showplot=T){
    fn <- warp_data$fn
    time <- warp_data$time
    qn <- warp_data$qn
    M <- nrow(qn)
    N <- ncol(qn)
    q0 <- warp_data$q0
    gam <- warp_data$gam
    # Set up for fPCA in q-space
    mq_new <- rowMeans(qn)
    m_new <- sign(fn[id,])*sqrt(abs(fn[id,]))  # scaled version
    mqn <- c(mq_new,mean(m_new))
    qn1 <- rbind(qn,m_new)

    # Calculate Vector Space of warping
    Tgam <- SqrtMean(gam)
    mu_psi <- Tgam$mu
    vec <- Tgam$vec

    # Joint fPCA --------------------------------------------------------------
    jointfPCAd <- function(qn, vec, C=1, m=3){
        M <- nrow(qn)
        N <- ncol(qn)
        time1 <- seq(0,2,length.out=M+nrow(vec))
        g <- rbind(qn,C*vec)
        mu_q <- rowMeans(qn)

        mu_g <- rowMeans(g)

        K <- cov(t(g))
        out.K <- svd(K, nu=m, nv=m)
        s <- out.K$d
        U <- out.K$u

        a <- matrix(0,N,m)
        for (i in 1:N){
            for (j in 1:m){
                a[i,j] <- (g[,i]-mu_g)%*%U[,j]
            }
        }

        qhat <- matrix(mu_q,M,N) + U[1:M,1:m] %*% t(a)
        vechat <- U[(M+1):length(time1),1:m] %*% t(a/C)
        psihat <- apply(vechat,2,function(x,mu_psi){exp_map(mu_psi,x)},mu_psi)
        gamhat <- apply(psihat,2,function(x,time){gam_mu = cumtrapz(time, x*x)
        gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))},
        seq(0,1,length.out=M-1))

        return(list(qhat=qhat,gamhat=gamhat,a=a,U=U[,1:m],s=s[1:m],mu_g=mu_g,cov=K,g=g))
    }


    # Find C ------------------------------------------------------------------
    findC <- function(C, qn, vec, q0, m){
        out.pca <- jointfPCAd(qn, vec, C, m)
        M <- nrow(qn)
        N <- ncol(qn)
        time <- seq(0, 1, length.out=M-1)

        d <- rep(0,N)
        for (i in 1:N){
            tmp <- warp_q_gamma(out.pca$qhat[1:(M-1),i], time, invertGamma(out.pca$gamhat[,i]))
            d[i] <- sum(trapz(time, (tmp-q0[,i])^2))
        }

        return(sum(d^2)/N)
    }

    m <- no
    if (is.null(C))
        C <- optimize(findC, c(0,1e4),qn=qn1,vec=vec,q0=q0,m=m)$minimum

    # Final PCA ---------------------------------------------------------------
    out.pca <- jointfPCAd(qn1, vec, C, m=m)

    # geodesic paths
    ci <- c(-1,0,1)
    q_pca <- array(0,dim=c(M,length(ci),m))
    f_pca <- array(0,dim=c(M,length(ci),m))
    N1 <- nrow(out.pca$U)
    for (j in 1:m){
        for (i in 1:length(ci)){
            qhat <- mqn + out.pca$U[1:(M+1),j] * ci[i]*sqrt(out.pca$s[j])
            vechat <- out.pca$U[(M+2):N1,j] * (ci[i]*sqrt(out.pca$s[j]))/C
            psihat <- exp_map(mu_psi,vechat)
            gamhat <- cumtrapz(seq(0,1,length.out=M), psihat*psihat)
            gamhat = (gamhat - min(gamhat))/(max(gamhat)-min(gamhat))
            if (sum(vechat)==0)
                gamhat <- seq(0,1,length.out = M)
            if (id == 1)
                fhat <- cumtrapz(time,qhat[1:M]*abs(qhat[1:M]))+sign(qhat[M+1])*(qhat[M+1]^2)
            else
                fhat <- cumtrapzmid(time,qhat[1:M]*abs(qhat[1:M]),sign(qhat[M+1])*(qhat[M+1]^2), id)
            f_pca[,i,j] <- warp_f_gamma(fhat, seq(0,1,length.out=M), gamhat)
            q_pca[,i,j] <- warp_q_gamma(qhat[1:M], seq(0,1,length.out=M), gamhat)
        }
    }

    jfpca <- list()
    jfpca$q_pca <- q_pca
    jfpca$f_pca <- f_pca
    jfpca$latent <- out.pca$s
    jfpca$coef <- out.pca$a
    jfpca$U <- out.pca$U
    jfpca$mu_psi <- mu_psi
    jfpca$mu_g <- out.pca$mu_g
    jfpca$id <- id
    jfpca$C <- C
    jfpca$time <- time
    jfpca$g <- out.pca$g
    jfpca$cov <- out.pca$cov

    class(jfpca) <- "jfpca"

    if (showplot){
        plot(jfpca)
    }

    return(jfpca)
}
