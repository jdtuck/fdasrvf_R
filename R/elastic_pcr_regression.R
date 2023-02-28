#' Elastic Linear Principal Component Regression
#'
#' This function identifies a regression model with phase-variability
#' using elastic pca
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param y vector of size \eqn{M} responses
#' @param time vector of size \eqn{N} describing the sample points
#' @param pca.method string specifying pca method (options = "combined",
#' "vert", or "horiz", default = "combined")
#' @param no scalar specify number of principal components (default=5)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel run in parallel (default = F)
#' @param C scale balance parameter for combined method (default = NULL)
#' @return Returns a pcr object containing
#' \item{alpha}{model intercept}
#' \item{b}{regressor vector}
#' \item{y}{response vector}
#' \item{warp_data}{fdawarp object of aligned data}
#' \item{pca}{pca object of principal components}
#' \item{SSE}{sum of squared errors}
#' \item{pca.method}{string specifying pca method used}
#' @keywords srvf alignment regression
#' @references J. D. Tucker, J. R. Lewis, and A. Srivastava, “Elastic
#'  Functional Principal Component Regression,” Statistical Analysis and Data
#'  Mining, 10.1002/sam.11399, 2018.
#' @export
elastic.pcr.regression <- function(f, y, time, pca.method="combined", no=5,
                                   smooth_data = FALSE, sparam = 25, parallel=F, C=NULL){

    pca.method1 <- pca.method
    pca.method <- pmatch(pca.method, c("combined", "vert", "horiz")) # 1 - combined, 2 - vert, 3 - horiz
    if (is.na(pca.method))
        stop("invalid method selection")

    if (smooth_data){
        f <- smooth.data(f,sparam)
    }

    N1 <- ncol(f)

    # Align Data --------------------------------------------------------------
    out <- time_warping(f, time, parallel = parallel)


    # Calculate PCA -----------------------------------------------------------
    if (pca.method == 1)
        out.pca <- jointFPCA(out, no, showplot=F, C=C)
    if (pca.method == 2)
        out.pca <- vertFPCA(out, no, showplot=F)
    if (pca.method == 3)
        out.pca <- horizFPCA(out ,no, showplot=F)

    # OLS using PCA basis
    lam <- 0  # regularization
    R <- 0
    Phi <- matrix(1, N1, no+1)
    Phi[,2:(no+1)] <- out.pca$coef
    xx <- t(Phi) %*% Phi
    inv_xx <- solve(xx + lam * R)
    xy <- t(Phi) %*% y
    b <- inv_xx %*% xy
    alpha <- b[1]
    b <- b[2:length(b)]

    # compute the SSE
    int_X <- rep(0, N1)
    for (ii in 1:N1){
        int_X[ii] <- sum(out.pca$coef[ii,] * b)
    }

    SSE <- sum((y-alpha-int_X)^2)

    out <- list(alpha=alpha, b=b, y=y, warp_data=out, pca=out.pca,
                SSE=SSE, pca.method=pca.method1)

    class(out) <- "pcr"

    return(out)
}
