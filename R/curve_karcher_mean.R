#' Karcher Mean of Curves
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic square-root velocity (srvf) framework.
#'
#' @param beta array (n,T,N) for N number of curves
#' @param mode Open ("O") or Closed ("C") curves
#' @param rotated Optimize over rotation (default = T)
#' @param maxit maximum number of iterations
#' @param ms string defining whether the Karcher mean ("mean") or Karcher median ("median") is returned (default = "mean")
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{betamean}{mean or median curve}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' \item{gam}{array of warping functions}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_karcher_mean(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
curve_karcher_mean <- function (beta, mode = "O", rotated = T, maxit = 20, ms = "mean") 
{
    if(ms!="mean"&ms!="median"){warning("ms must be either \"mean\" or \"median\". ms has been set to \"mean\"",immediate. = T)}
    if(ms!="median"){ms = "mean"}
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]
    q = array(0, c(n, T1, N))
    for (ii in 1:N) {
        q[, , ii] = curve_to_q(beta[, , ii])
    }
    mnq = rowMeans(q[1, , ])
    dqq = sqrt(colSums((q[1, , ] - matrix(mnq, ncol = N, nrow = T1))^2))
    min_ind = which.min(dqq)
    mu = q[, , min_ind]
    betamean = beta[, , min_ind]
    delta = 0.5
    tolv = 1e-04
    told = 5 * 0.001
    itr = 1
    sumd = rep(0, maxit + 1)
    v = array(0, c(n, T1, N))
    normvbar = rep(0, maxit + 1)
    if(ms == "median"){ #run for median only, saves memory if getting mean
        d_i = rep(0,N) #include vector for norm calculations
        v_d = array(0, c(n, T1, N)) #include array to hold v_i / d_i
    }
    while (itr < maxit) {
        cat(sprintf("Iteration: %d\n", itr))
        mu = mu/sqrt(innerprod_q2(mu, mu))
        
        for (i in 1:N) {
            out = karcher_calc(beta[, , i], q[, , i], betamean, 
                               mu, rotated, mode)
            v[, , i] = out$v
            if(ms == "median"){ #run for median only, saves computation time if getting mean
                d_i[i] = sqrt(innerprod_q2(v[,,i], v[,,i])) #calculate sqrt of norm of v_i
                v_d[,,i] = v[,,i]/d_i[i] #normalize v_i
            }
            sumd[itr + 1] = sumd[itr + 1] + out$d^2
        }
        if(ms == "median"){#run for median only
            sumv = rowSums(v_d, dims = 2)
            sum_dinv = sum(1/d_i)
            vbar = sumv/sum_dinv
        }
        else{ #run for mean only
            sumv = rowSums(v, dims = 2)
            vbar = sumv/N
        }
        normvbar[itr] = sqrt(innerprod_q2(vbar, vbar))
        normv = normvbar[itr] 
        if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > 
                               told)) {
            mu = cos(delta * normvbar[itr]) * mu + sin(delta * 
                                                           normvbar[itr]) * vbar/normvbar[itr]
            if (mode == "C") {
                mu = project_curve(mu)
            }
            x = q_to_curve(mu)
            a = -1 * calculatecentroid(x)
            dim(a) = c(length(a), 1)
            betamean = x + repmat(a, 1, T1)
        }
        else {
            break
        }
        itr = itr + 1
    }
    ifelse(ms=="median",type<-"Karcher Median",type<-"Karcher Mean")
    return(list(mu = mu, type = type, betamean = betamean, v = v, q = q))
}
