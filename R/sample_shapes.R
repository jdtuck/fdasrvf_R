#' Sample shapes from model
#'
#' @param x An object of class `fdacurve` typically produced by [multivariate_karcher_mean()]
#' @param no number of principal components
#' @param numSamp number of samples
#' @return Returns a `fdacurve` object containing \item{betans}{random aligned curves}
#' \item{qns}{random aligned srvfs}
#' \item{betans}{random curves}
#' \item{qs}{random srvfs}
#' \item{gams}{random reparameterization functions}
#' \item{R}{random rotation matrices}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'    Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- multivariate_karcher_mean(beta[, , 1, 1:5], scale=TRUE, maxit = 2)
#' # note: use more shapes, small for speed
#' out.samples <- sample_shapes(out)
sample_shapes <- function(x, no=3, numSamp=10){

    mode = x$mode

    if (!x$scale)
      cli::cli_abort("Not implemented for `scale = FALSE`. Please set `scale = TRUE` in multivariate_karcher_mean.")

    K <- multivariate_karcher_cov(x)
    mu <- x$mu
    mu <- mu/sqrt(innerprod_q2(mu, mu))

    n = nrow(mu)
    T1 = ncol(mu)

    out = svd(K)
    U = out$u
    s = out$d
    V = out$v

    if (mode == "O"){
        N = 2
    } else {
        N = 5  #10
    }

    epsilon = 1./(N-1)

    q1 = mu
    q2 = mu
    samples = array(0,dim=c(n,T1,numSamp))
    samples.q = array(0,dim=c(n,T1,numSamp))

    # distribution if scales
    scale = rep(1, numSamp)
    scale_min = min(x$len_q)
    scale_max = max(x$len_q)
    if (!x$scale){
      scale = stats::runif(numSamp,scale_min,scale_max)
    }

    for (i in 1:numSamp){
        v = matrix(0, 2, T1)
        for (m in 1:no){
            v = v + stats::rnorm(1)*sqrt(s[m])*c(U[1:T1,m], U[(T1+1):(2*T1),m])
        }

        q1 = mu
        for (j in 1:(N-1)){
            normv = sqrt(innerprod_q2(v,v))

            if (normv < 1e-4){
                q2 = mu
            } else {
                q2 = cos(epsilon*normv)*q1+sin(epsilon*normv)*v/normv
                if (mode == "C"){
                    q2 = project_curve(q2)
                }
            }

            # Parallel translate tangent vector
            # basis2 = find_basis_normal(q2)
            # v = parallel_translate(v, q1, q2, basis2, mode)

            # q1 = q2
        }

        beta = q_to_curve(q2, scale[i])
        centroid = calculatecentroid(beta)
        dim(centroid) = c(length(centroid),1)
        beta = beta - repmat(centroid,1,T1)
        samples[,,i] = beta
        samples.q[,,i] = curve_to_q(beta, x$scale)$q
    }

    # random warping generation
    rgam = randomGamma(x$gam,numSamp)
    gams = matrix(0,numSamp,T1)
    for (k in 1:numSamp){
      gams[k,] = invertGamma(rgam[k,])
    }
    gams = t(gams)

    # random rotation angles
    if (x$rotation){
      # currently for R^2
      N = dim(x$R)[3]
      theta = rep(0, N)
      for (i in 1:N){
        theta[i] = acos(x$R[1,1,i])
      }
      mu_theta = mean(theta)
      sd_theta = stats::sd(theta)
      R = array(0, dim=c(2,2,numSamp))
      for (k in 1:numSamp){
        theta = sample_vonmises(mu_theta, sd_theta)
        R[,,k] = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)
      }
    } else {
      R = array(0, dim=c(2,2,numSamp))
      for (k in 1:numSamp){
        theta = 0
        R[,,k] = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),2,2)
      }
    }

    samples1 = samples
    samples.q1 = samples.q
    for (k in 1:numSamp){
      beta1 = R[,,k]%*%samples1[,,k]
      samples1[,,k] = group_action_by_gamma_coord(beta1, gams[,k])
      samples.q1[,,k] = curve_to_q(samples1[,,k], x$scale)$q

    }

    x$rsamps = TRUE
    x$betans = samples
    x$betas = samples1
    x$qns = samples.q
    x$qs = samples.q1
    x$gams = gams
    x$R = R

    return(x)
}


sample_vonmises <- function(mu, sigma){
  while (TRUE){
    u1 = stats::runif(1)
    u2 = stats::runif(1)
    phi = 2 * pi * u1

    R = exp((cos(phi-mu) - 1)/ sigma^2)

    if (u2 < R){
      break
    }
  }

    return(phi)
}
