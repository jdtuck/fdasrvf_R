#' Karcher Mean of Curves
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic square-root velocity (srvf) framework.
#'
#' @param beta array (n,T,N) for N number of curves
#' @param mode Open ("O") or Closed ("C") curves
#' @param rotated Optimize over rotation (default = T)
#' @param scale Include scale (default = F)
#' @param maxit maximum number of iterations
#' @param ms string defining whether the Karcher mean ("mean") or Karcher median ("median") is returned (default = "mean")
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{beta}{centered data}
#' \item{betamean}{mean or median curve}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' \item{gam}{array of warping functions}
#' \item{cent}{centers of original curves}
#' \item{len}{length of curves}
#' \item{len_q}{length of srvfs}
#' \item{mean_scale}{mean length}
#' \item{mean_scale_q}{mean length srvf}
#' \item{E}{energy}
#' \item{qun}{cost function}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more shapes, small for speed
curve_karcher_mean <- function (beta, mode = "O", rotated = T, scale = F, maxit = 20, ms = "mean")
{
    if(ms!="mean"&ms!="median"){warning("ms must be either \"mean\" or \"median\". ms has been set to \"mean\"",immediate. = T)}
    if(ms!="median"){ms = "mean"}

    mean_scale=NA
    mean_scale_q=NA
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]
    q = array(0, c(n, T1, N))
    len = rep(0,N)
    len_q = rep(0,N)
    cent = matrix(0,n,N)
    for (ii in 1:N) {
        beta1 = beta[,,ii]
        centroid1 = calculatecentroid(beta1)
        cent[,ii] = -1*centroid1
        dim(centroid1) = c(length(centroid1),1)
        beta1 = beta1 - repmat(centroid1,1,T1)
        beta[,,ii] = beta1
        out = curve_to_q(beta1)
        q[, , ii] = out$q
        len[ii] = out$len
        len_q[ii] = out$lenq
    }

    mu = q[, , 1]
    bmu = beta[, , 1]
    delta = 0.5
    tolv = 1e-04
    told = 5 * 0.001
    itr = 1
    sumd = rep(0, maxit + 1)
    sumd[1] = Inf
    v = array(0, c(n, T1, N))
    normvbar = rep(0, maxit + 1)
    if(ms == "median"){ #run for median only, saves memory if getting mean
        d_i = rep(0,N) #include vector for norm calculations
        v_d = array(0, c(n, T1, N)) #include array to hold v_i / d_i
    }


    cat("\nInitializing...\n")
    gam = matrix(0,T1,N)
    for (k in 1:N) {
        out = find_rotation_seed_unqiue(mu,q[, , k],mode)
        gam[,k] = out$gambest
    }

    gam = t(gam)
    gamI = SqrtMeanInverse(t(gam))
    bmu = group_action_by_gamma_coord(bmu, gamI)
    mu = curve_to_q(bmu)$q
    mu[is.nan(mu)] <- 0

    while (itr < maxit) {
        cat(sprintf("Iteration: %d\n", itr))
        mu = mu/sqrt(innerprod_q2(mu, mu))

        if (mode=="C"){
            basis = find_basis_normal(mu)
        }

        for (i in 1:N) {
            q1 = q[, , i]

            out = find_rotation_seed_unqiue(mu,q1,mode)
            qn_t = out$q2best/sqrt(innerprod_q2(out$q2best,out$q2best))

            q1dotq2 = innerprod_q2(mu,qn_t)

            if (q1dotq2 > 1){
                q1dotq2 = 1
            }
            if (q1dotq2 < -1){
                q1dotq2 = -1
            }

            dist = acos(q1dotq2)

            u = qn_t - q1dotq2 * q1
            normu = sqrt(innerprod_q2(u,u))
            if (normu > 1e-4){
                w = u*acos(q1dotq2)/normu
            } else {
                w = matrix(0, nrow(beta1), T1)
            }

            if (mode=="O"){
                v[, , i] = w
            } else {
                v[, , i] = project_tangent(w, q1, basis)
            }

            if(ms == "median"){ #run for median only, saves computation time if getting mean
                d_i[i] = sqrt(innerprod_q2(v[,,i], v[,,i])) #calculate sqrt of norm of v_i
                if (d_i[i]>0){
                    v_d[,,i] = v[,,i]/d_i[i] #normalize v_i
                } else{
                    v_d[,,i] = v[,,i]
                }

            }
            sumd[itr + 1] = sumd[itr + 1] + dist^2
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
        if ((sumd[itr]-sumd[itr+1]) < 0){
            break
        } else if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > told)) {
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

    if (scale){
        mean_scale = prod(len)^(1/length(len))
        mean_scale_q = prod(len_q)^(1/length(len))
        betamean = mean_scale*betamean
    }

    ifelse(ms=="median",type<-"Karcher Median",type<-"Karcher Mean")
    return(list(beta = beta, mu = mu, type = type, betamean = betamean, v = v, q = q,
                E=normvbar[1:itr], cent = cent, len = len, len_q = len_q,
                qun = sumd[1:itr], mean_scale = mean_scale, mean_scale_q=mean_scale_q))
}
