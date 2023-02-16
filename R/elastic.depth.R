#' Calculates elastic depth
#'
#' This functions calculates the elastic depth between set of functions
#'
#' @param f matrix of N function of M time points (MxN)
#' @param time sample points of functions
#' @param lambda controls amount of warping (default = 0)
#' @param parallel run computation in parallel (default = T)
#' @return Returns a list containing \item{amp}{amplitude depth}
#' \item{phase}{phase depth}
#' @keywords depth
#' @concept srvf alignment
#' @references T. Harris, J. D. Tucker, B. Li, and L. Shand, "Elastic depths for detecting shape anomalies in functional data," Technometrics, 10.1080/00401706.2020.1811156, 2020.
#' @export
#' @examples
#' depths <- elastic.depth(simu_data$f[, 1:4], simu_data$time)
elastic.depth <- function(f,time,lambda = 0, parallel = FALSE){
    if (parallel){
        cores = detectCores()-1
        cl = makeCluster(cores)
        registerDoParallel(cl)
    } else
    {
        registerDoSEQ()
    }

    obs = nrow(f)
    fns = ncol(f)

    amp_dist = matrix(0, fns, fns)
    phs_dist = matrix(0, fns, fns)
    k = 0

    for (f1 in 1:(fns-1)) {

        dist<-foreach(k = f1:ncol(f), .combine=cbind,.packages='fdasrvf') %dopar% {
            out = elastic.distance(f[,f1], f[,k], time)

            list(out$Dy,out$Dx)
        }

        N = ncol(f)-f1+1
        phs = unlist(dist[2,])
        dim(phs)=c(1,N)
        amp = unlist(dist[1,])
        dim(amp)=c(1,N)

        phs_dist[f1, f1:fns] = phs
        amp_dist[f1, f1:fns] = amp
    }

    amp_dist = amp_dist + t(amp_dist)
    phs_dist = phs_dist + t(phs_dist)

    amp = 1 / (1 + apply(amp_dist, 1, median))
    phase = 1 / (1 + apply(phs_dist, 1, median))
    phase = ((2+pi)/pi) * (phase - 2/(2+pi))

    if (parallel){
        stopCluster(cl)
    }

    return(list(amp=amp,phase=phase))
}
