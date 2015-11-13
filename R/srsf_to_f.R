#' Convert SRVF tp f
#'
#' This function converts srvfs to functions
#'
#' @param q matrix of srvfs
#' @param time time
#' @param f0: initial value of f
#' @return f matrix of functions
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(simu_data$f,simu_data$time)
#' f = f_to_srvf(q,simu_data$time,f[1,])
srsf_to_f <- function(q,time,f0=0.0){
    if is.null(dim(q)){
        integrand = q*abs(q)
        f = f0 + cumtrapz(time, integrand)
    } else {
        M = nrows(q);
        N = ncols(q);
        f = matrix(0,M,N)
        for (i in 1:N){
            qnorm = abs(q[,i])
            integrand = q[,i]*qnorm
            f[,i] = f0[i] + cumtrapz(time, integrand)
        }
    }

    return(f)
}
