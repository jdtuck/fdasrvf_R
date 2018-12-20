#' Elastic Logisitc Prinipcal Component Regression
#'
#' This function identifies a logistic regression model with phase-variability
#' using elastic pca
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param y vector of size \eqn{M} lables
#' @param time vector of size \eqn{N} describing the sample points
#' @param pca.method string specifing pca method (options = "combined",
#' "vert", or "horiz", default = "combined")
#' @param no scalar specifify number of principal components (default=5)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @return Returns a lpcr object containing
#' \item{alpha}{model intercept}
#' \item{b}{regressor vector}
#' \item{y}{label vector}
#' \item{warp_data}{fdawarp object of aligned data}
#' \item{pca}{pca object of principal components}
#' \item{Loss}{logistic loss}
#' \item{pca.method}{string specifing pca method used}
#' @keywords srvf alignment regression
#' @references J. D. Tucker, J. R. Lewis, and A. Srivastava, “Elastic
#'  Functional Principal Component Regression,” Statistical Analysis and Data
#'  Mining, 10.1002/sam.11399, 2018.
#' @export
elastic.lpcr.regression <- function(f, y, time, pca.method="combined", no=5,
                                    smooth_data = FALSE, sparam = 25){

    pca.method1 <- pca.method
    pca.method <- pmatch(pca.method, c("combined", "vert", "horiz")) # 1 - combined, 2 - vert, 3 - horiz
    if (is.na(pca.method))
        stop("invalid method selection")

    if (smooth_data){
        f = smooth.data(f,sparam)
    }

    N1 <- ncol(f)

    # Align Data --------------------------------------------------------------
    out <- time_warping(f, time, parallel = T, showplot = F)


    # Calculate PCA -----------------------------------------------------------
    if (pca.method == 1)
        out.pca <- jointFPCA(out, no, showplot=F)
    if (pca.method == 2)
        out.pca <- vertFPCA(out, no, showplot=F)
    if (pca.method == 3)
        out.pca <- horizFPCA(out ,no, showplot=F)

    # OLS using PCA basis
    lam <- 0  # regularization
    R <- 0
    Phi <- matrix(1, N1, no+1)
    Phi[,2:(no+1)] <- out.pca$coef
    # Find alpha and beta using l_bfgs
    b0 <- rep(0,no+1)
    out1 <- optim(b0, logit_loss, gr = logit_gradient, Phi, y,
                  method = "L-BFGS-B", control = list(maxit=200,pgtol=1e-10))
    b <- out1$par

    # compute the Loss
    LL = logit_loss(b,Phi,y)

    alpha<-b[1]
    b <- b[2:length(b)]

    out <- list(alpha=alpha, b=b, y=y, warp_data=out, pca=out.pca,
                LL=LL, pca.method=pca.method1)

    class(out) <- "lpcr"

    return(out)
}
