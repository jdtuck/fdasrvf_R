#' Horizontal Functional Principal Component Analysis
#'
#' This function calculates vertical functional principal component analysis
#' on aligned data
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param no number of principal components to extract
#' @param var_exp compute no based on value percent variance explained (example: 0.95)
#'                will override `no`
#' @param ci geodesic standard deviations (default = c(-1,0,1))
#' @param showplot show plots of principal directions (default = T)
#' @return Returns a hfpca object containing \item{gam_pca}{warping functions principal directions}
#' \item{psi_pca}{srvf principal directions}
#' \item{latent}{latent values}
#' \item{U}{eigenvectors}
#' \item{vec}{shooting vectors}
#' \item{mu}{Karcher Mean}
#' @keywords srvf alignment
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' hfpca <- horizFPCA(simu_warp, no = 3)
horizFPCA <- function(warp_data, no=3, var_exp=NULL, ci=c(-1,0,1),
                      showplot = TRUE){
    gam <- warp_data$warping_functions
    tmp = SqrtMean(gam)
    vec = tmp$vec
    mu = tmp$mu

    # TFPCA
    K = stats::cov(t(vec)) #out$sigma

    out = svd(K)
    s = out$d
    U = out$u
    TT = nrow(vec) + 1
    vm = rowMeans(vec)

    # Parameters
    if (!is.null(var_exp)){
      cumm_coef <- cumsum(s)/sum(s)
      tmp = which(cumm_coef <= var_exp)
      no = tmp[length(tmp)]
    }

    no_pca = 1:no

    gam_pca = array(0,dim=c(length(ci),length(mu)+1,no))
    psi_pca = array(0,dim=c(length(ci),length(mu),no))
    for (j in no_pca){
        cnt = 1
        for (k in ci){
            v = k*sqrt(s[j])*U[,j]
            vn = pvecnorm(v,2)/sqrt(TT)
            if (vn < 0.0001){
                psi_pca[cnt,,j] = mu
            }else{
                psi_pca[cnt,,j] = cos(vn)*mu+sin(vn)*v/vn
            }
            tmp = rep(0,TT)
            tmp[2:TT] = cumsum(psi_pca[cnt,,j]*psi_pca[cnt,,j])
            gam_pca[cnt,,j] = (tmp-tmp[1]) / (tmp[length(tmp)] - tmp[1])
            cnt = cnt + 1
        }
    }

    N2 = dim(gam)[2]
    c = matrix(0,N2,no)
    for (k in no_pca){
      for (i in 1:N2){
        c[i,k] = sum((vec[,i]-vm)*U[,k])
      }
    }

    hfpca = list()
    hfpca$gam_pca = gam_pca
    hfpca$psi_pca = psi_pca
    hfpca$latent = s[no_pca]
    hfpca$U = U[,no_pca]
    hfpca$coef = c[,no_pca]
    hfpca$vec = vec
    hfpca$mu = mu
    hfpca$vm = vm
    hfpca$warp_data = warp_data

    class(hfpca) <- "hfpca"

    if (showplot){
        plot(hfpca)
    }
    return(hfpca)
}
