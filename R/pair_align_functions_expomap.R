

# TODO:
# objects don't need to pass around lists with time (x), only the function values (y); in fact, can get rid of most of time references
# don't order all the time! sort f1, f2 once at the outset
# integrate with fdasrvf functions
# spacing & formatting

################################################################################
# calculate exp1(g), exp1inv(psi)
# g, psi: function in the form of list$x, list$y
# returns exp1(g) or exp1inv(psi), function in the form of list$x, list$y
################################################################################

# new function
f.exp1 <- function(g) {
  return(list(x = g$x,
              y = calcY(f.L2norm(g), g$y)))
}

f.exp1inv <- function(psi) {
  #psi=onef
  #psi=f.phi(test)
  x = psi$x
  y = psi$y

  #intergral estimation
  ind = order(x)
  inner = round(trapzCpp(x[ind], y[ind]), 10)

  if ((inner < 1.001) && (inner >= 1)) {
    inner = 1
  }
  if ((inner <= -1) && (inner > -1.001)) {
    inner = -1
  }
  if ((inner < (-1)) ||
      (inner > 1)) {
    print(paste("exp1inv: can't calculate the acos of:", inner))
  }

  theta = acos(inner)
  yy = theta / sin(theta) * (y - rep(inner, length(y)))

  if (theta == 0) {
    yy = rep(0, length(x))
  }

  return(list(x = x, y = yy))
}

# function for calculating the next MCMC sample given current state
f.updateg.pw <-
  function(g.coef.curr,
           g.basis,
           var1.curr,
           q1,
           q2,
           SSE.curr = 0,
           propose.g.coef) {
    g.coef.prop = propose.g.coef(g.coef.curr = g.coef.curr)
    while (min(f.exp1(
      f.basistofunction(
        f.domain = g.basis$x,
        coef = g.coef.prop$prop,
        basis = g.basis
      )
    )$y) < 0) {
      g.coef.prop = propose.g.coef(g.coef.curr = g.coef.curr)
    }

    if (SSE.curr == 0) {
      SSE.curr = f.SSEg.pw(
        g = f.basistofunction(
          f.domain = g.basis$x,
          coef = g.coef.curr,
          basis = g.basis
        ),
        q1 = q1,
        q2 = q2
      )
    }

    SSE.prop = f.SSEg.pw(
      g = f.basistofunction(
        f.domain = g.basis$x,
        coef = g.coef.prop$prop,
        basis = g.basis
      ),
      q1 = q1,
      q2 = q2
    )

    logl.curr = f.logl.pw(
      g = f.basistofunction(
        f.domain = g.basis$x,
        coef = g.coef.curr,
        basis = g.basis
      ),
      var1 = var1.curr,
      q1 = q1,
      q2 = q2,
      SSEg = SSE.curr
    )
    logl.prop = f.logl.pw(
      g = f.basistofunction(
        f.domain = g.basis$x,
        coef = g.coef.prop$prop,
        basis = g.basis
      ),
      var1 = var1.curr,
      q1 = q1,
      q2 = q2,
      SSEg = SSE.prop
    )

    ratio = min(1, exp(logl.prop - logl.curr))  #prob to accept the proposal
    u = runif(1)
    if (u <= ratio) {
      return(
        list(
          g.coef = g.coef.prop$prop,
          logl = logl.prop,
          SSE = SSE.prop,
          accept = TRUE,
          zpcnInd = g.coef.prop$ind
        )
      )
    }
    if (u > ratio) {
      return(
        list(
          g.coef = g.coef.curr,
          logl = logl.curr,
          SSE = SSE.curr,
          accept = FALSE,
          zpcnInd = g.coef.prop$ind
        )
      )
    }
  }

################################################################################
# For pairwise registration, evaluate the loglikelihood of g, given q1 and q2
# g, q1, q2 are all given in the form of list$x and list$y
# var1: model variance
# SSEg: if not provided, than re-calculate
# returns a numeric value which is logl(g|q1,q2), see Eq 10 of JCGS
################################################################################
#SSEg: sum of sq error= sum over ti of { q1(ti)-{q2,g}(ti) }^2 (Eq 11 of JCGS)

f.SSEg.pw <- function(g, q1, q2) {
  par.domain = g$x
  obs.domain = q1$x
  exp1g.temp = f.predictfunction(f = f.exp1(g), at = obs.domain)
  pt <- c(0, cuL2norm2(x = obs.domain, y = exp1g.temp$y))
  vec = (q1$y - f.predictfunction(f = q2, at = pt)$y * (exp1g.temp$y)) ^
    2
  return(sum(vec))
}

f.logl.pw <- function(g, q1, q2, var1, SSEg = 0) {
  if (SSEg == 0) {
    SSEg = f.SSEg.pw(g = g, q1 = q1, q2 = q2)
  }
  n = length(q1$y)
  return(n * log(1 / sqrt(2 * pi)) - n * log(sqrt(var1)) - SSEg / (2 * var1))
}



################################################################################
# calculate Q(f), Qinv(q)
# f,q: function in the form of list$x, list$y
# fini: f(0)
# returns Q(f),Qinv(q), function in the form of list$x, list$y
################################################################################
f.Q <- function(f) {
  d = f.predictfunction(f, at = f$x, deriv = 1)$y
  result = list(x = f$x, y = sign(d) * sqrt(abs(d)))
  return(result)
}

f.Qinv <- function(q, fini = 0) {
  result = rep(NA, length(q$x))
  for (i in 1:length(result)) {
    y = q$y[1:i]
    x = q$x[1:i]
    ind = order(x)
    result[i] = trapzCpp(x[ind], (y[ind] * abs(y[ind])))
  }
  result = result + fini
  result[1] = fini
  return(list(x = q$x, y = result))
}



################################################################################
# Extrapolate a function given by a discrete vector
# f: function in the form of list$x, list$y
# at: t values for which f(t) is returned
# deriv: can calculate derivative
# method: smoothing method: 'linear' (default), 'cubic'
# returns: $y==f(at), $x==at
################################################################################
f.predictfunction <- function(f,
                              at,
                              deriv = 0,
                              method = "linear") {
  if (method == "cubic") {
    result = predict(smooth.spline(
      f$x,
      f$y,
      all.knots = F,
      nknots = (length(f$x) - 10)
    ), at, deriv = deriv)
    if ((at[1] == f$x[1]) && (deriv == 0)) {
      result$y[1] = f$y[1]
    }
  }
  if (method == "linear") {
    if (deriv == 0) {
      result = approx(f$x, f$y, xout = at, rule = 2)
    }
    if (deriv == 1) {
      #f=list(x=c(0,0.5,1),y=c(1,2,4))
      #at=c(0,0.25,0.5,0.75,1)
      fmod = approx(f$x, f$y, rule = 2, xout = at)
      diffy1 = c(0, diff(fmod$y))
      diffy2 = c(diff(fmod$y), 0)
      diffx1 = c(0, diff(fmod$x))
      diffx2 = c(diff(fmod$x), 0)

      (diffy2 + diffy1) / (diffx2 + diffx1)
      result = list(x = at,
                    y = (diffy2 + diffy1) / (diffx2 + diffx1))

    }
    if (deriv >= 2) {
      stop("f.predictfunction: can't calculate >=2nd derivative using linear method")
    }
  }
  return (result)
}



################################################################################
# calculate L2 norm of a function, using trapezoid rule for integration
# f:function in the form of list$x, list$y
# returns ||f||, a numeric value
################################################################################

f.L2norm <- function(f) {
  return(order_l2norm(f$x, f$y))
}


################################################################################
# Different basis functions b_i()
# f.domain: grid on which b_i() is to be returned
# numBasis: numeric value, number of basis functions used (note: #basis = #coef/2 for Fourier basis)
#fourier.p: period of the Fourier basis used
# returns a list:
#     $matrix: with nrow=length(t) and ncol=numBasis (or numBasis*2 for Fourier)
#     $x: f.domain
################################################################################
basis.fourier <- function(f.domain = 0,
                          numBasis = 1,
                          fourier.p = 1) {
  result = matrix(NA, nrow = length(f.domain), ncol = 2 * numBasis)
  for (i in 1:(2 * numBasis)) {
    j = ceiling(i / 2)
    if ((i %% 2) == 1) {
      result[, i] = sqrt(2) * sin(2 * j * pi * f.domain / fourier.p)
    }
    if ((i %% 2) == 0) {
      result[, i] = sqrt(2) * cos(2 * j * pi * f.domain / fourier.p)
    }
  }
  return(list(x = f.domain, matrix = result))
}

################################################################################
# Given the coefficients of basis functions, returns the actual function on a grid
# f.domain: numeric vector, grid of the actual function to return
# coefconst: leading constant term
# coef: numeric vector, coefficients of the basis functions
#       Note: if #coef < #basis functions, only the first #coef basis functions will be used
# basis: in the form of list$x, list$matrix
# plotf: if true, show a plot of the function generated
# returns the generated function in the form of list$x=f.domain, list$y
################################################################################
f.basistofunction <-
  function(f.domain,
           coefconst = 0,
           coef,
           basis,
           plotf = F) {
    if (dim(basis$matrix)[2] < length(coef)) {
      stop('In f.basistofunction, #coeffients exceeds #basis functions.')
    }
    result = list(x = basis$x,
                  y = basis$matrix[, (1:length(coef))] %*% matrix(coef) + coefconst)
    result = f.predictfunction(f = result, at = f.domain)
    if (plotf) {
      plot(result)
    }
    return(result)
  }


############################################################################################
# Calculate exp_psi(g), expinv_psi(psi2)
# g, psi: function in the form of list$x, list$y
# returns exp_psi(g) or expinv_psi(psi2), function in the form of list$x, list$y
############################################################################################
f.exppsi <- function(psi, g) {
  #psi=list(x=f.domain,y=MCMC.sim.psi$MCMC[,1])
  #g=list(x=f.domain,y=MCMC.sim$MCMC[,1])
  area = round(f.L2norm(g), 10)
  y = cos(area) * (psi$y) + sin(area) / area * (g$y)
  if (area == 0) {
    y = psi$y
  }
  return(list(x = (g$x), y = y))
}

f.exppsiinv <- function(psi, psi2) {
  #psi=list(x=f.domain, y=rep(1,length(f.domain)))
  #psi2=list(x=f.domain,y=MCMC.sim.psi$MCMC[,2])
  x = psi$x

  #intergral estimation
  ind = order(x)
  inner = round(trapz(x[ind], (psi$y)[ind] * (psi2$y)[ind]), 10)
  if ((inner < 1.001) && (inner >= 1)) {
    inner = 1
  }
  if ((inner < 1.05) && (inner >= 1.001)) {
    print(paste("exppsiinv: caution! acos of:", inner, "is set to 1..."))
    inner = 1
  }
  if ((inner <= -1) && (inner > -1.001)) {
    inner = -1
  }
  if ((inner <= -1.001) && (inner > -1.05)) {
    print(paste("exppsiinv: caution! acos of:", inner, "is set to -1..."))
    inner = -1
  }
  if ((inner < (-1)) ||
      (inner > 1)) {
    print(paste("exppsiinv: can't calculate the acos of:", inner))
  }
  theta = acos(inner)
  yy = theta / sin(theta) * ((psi2$y) - inner * (psi$y))
  if (theta == 0) {
    yy = rep(0, length(x))
  }
  return(list(x = x, y = yy))
}


############################################################################################
# Calculate Karcher mean/median with Alg1/Alg2 in (Kurtek,2014)
# x: vector of length = length(domain of the psi's)
# y: M columns, each of length = length(x)
# e1, e2: small positive constants
# method: 'ext' = extrinsic, 'int' = intrinsic
# returns posterier mean/median of M psi's (of form $x, $y)
############################################################################################
f.psimean <- function(x,
                      y,
                      e1 = 0.001,
                      e2 = 0.3,
                      method = 'ext') {
  if (method == 'int') {
    p = dim(y)[2]
    vM = matrix(NA, nrow = length(x), ncol = p)
    #energy=numeric()

    #initialize
    result = rowMeans(y[, ]) / f.L2norm(list(x = x, y = rowMeans(y[, ])))
    vbar = rep(1, length(x))
    j = 1
    while (f.L2norm(list(x = x, y = vbar)) >= e1) {
      for (i in 1:p) {
        #i=1
        vM[, i] = f.exppsiinv(psi = list(x = x, y = result),
                              psi2 = list(x = x, y = y[, i]))$y
        #dM[i]=acos(trapz(f.domain, result*psiM[,i]))
      }
      #energy=c(energy,sum(dM^2))
      vbar = rowMeans(vM)
      #update
      result = f.exppsi(psi = list(x = x, y = result),
                        g = list(x = x, y = e2 * vbar))$y
      j = j + 1
      #if((j%%5)==0){print(paste("updated step",j,"... Norm is", f.L2norm(list(x=x,y=vbar)) ))}
    }
  }

  if (method == 'ext') {
    rmy <- rowMeans(y)
    result <- rmy / f.L2norm(list(x = x, y = rmy))
  }

  return(list(x = x, y = result))
}


############################################################################################
# calculate phi(gamma), phiinv(psi)
# gamma, psi: function in the form of list$x, list$y
# returns phi(gamma) or phiinv(psi), function in the form of list$x, list$y
############################################################################################
f.phi <- function(gamma) {
  f.domain = gamma$x
  k = f.predictfunction(gamma, at = f.domain, deriv = 1)$y
  if (length(which(k < 0)) != 0) {
    #print(i)

    #print(paste("phi: derivative of gamma at",f.domain[which(k<0)],"is",k[which(k<0)],"; set to 0"))
    #plot(gamma$x,gamma$y)
    #abline(v=f.domain[which(k<0)],lty='dotted')
    k[which(k < 0)] = 0
  }
  result = list(x = f.domain, y = sqrt(k))
  if (f.L2norm(result) >= (1.01) || f.L2norm(result) <= (0.99)) {
    #print(paste("phi: L2norm of phi(gamma) is",f.L2norm(result),"; set to 1"))
  }
  result$y = (result$y) / f.L2norm(result)
  return(result)
}
## the function returned by phi(gamma) = psi is always positive and has L2norm 1

f.phiinv <- function(psi) {
  f.domain = psi$x

  result <- c(0, cuL2norm2(x = f.domain, y = psi$y))

  #  result=rep(NA,length(f.domain))
  #  result[1]=0
  #  for(i in 2:length(result)){
  #    result[i]=f.L2norm( list(x=f.domain[1:i],y=psi$y[1:i]) )^2
  #  }

  return(list(x = f.domain, y = result))
}




#' Align two functions using geometric properties of warping functions
#'
#' This function aligns two functions using Bayesian framework. It will align
#' f2 to f1. It is based on mapping warping functions to a hypersphere, and a
#' subsequent exponential mapping to a tangent space. In the tangent space,
#' the Z-mixture pCN algorithm is used to explore both local and global
#' structure in the posterior distribution.
#'
#' The Z-mixture pCN algorithm uses a mixture distribution for the proposal
#' distribution, controlled by input parameter zpcn. The zpcn$betas must be
#' between 0 and 1, and are the coefficients of the mixture components, with
#' larger coefficients corresponding to larger shifts in parameter space. The
#' zpcn$probs give the probability of each shift size.
#'
#' @param f1 observed data, numeric vector
#' @param f2 observed data, numeric vector
#' @param timet sample points of functions
#' @param iter length of the chain
#' @param burnin number of burnin MCMC iterations
#' @param alpha0,beta0 IG parameters for the prior of sigma1
#' @param zpcn list of mixture coefficients and prior probabilities for
#' Z-mixture pCN algorithm of the form list(betas, probs), where betas and
#' probs are numeric vectors of equal length
#' @param propvar variance of proposal distribution
#' @param init.coef initial coefficients of warping function in exponential map;
#' length must be even
#' @param npoints number of sample points to use during alignment
#' @param extrainfo T/F whether additional information is returned
#' @return Returns a list containing
#' \item{f2_warped}{f2 aligned to f1}
#' \item{gamma}{Posterior mean gamma function}
#' \item{g.coef}{matrix with iter columns, posterior draws of g.coef}
#' \item{psi}{Posterior mean psi function}
#' \item{sigma1}{numeric vector of length iter, posterior draws of sigma1}
#' \item{accept}{Boolean acceptance for each sample (if extrainfo=TRUE)}
#' \item{betas.ind}{Index of zpcn mixture component for each sample (if extrainfo=TRUE)}
#' \item{logl}{numeric vector of length iter, posterior loglikelihood (if extrainfo=TRUE)}
# #' \item{psi.mat}{Matrix of all posterior draws of psi (if extrainfo=TRUE)}
#' \item{gamma_mat}{Matrix of all posterior draws of gamma (if extrainfo=TRUE)}
#' \item{gamma_q025}{Lower 0.025 quantile of gamma (if extrainfo=TRUE)}
#' \item{gamma_q975}{Upper 0.975 quantile of gamma (if extrainfo=TRUE)}
#' \item{sigma_eff_size}{Effective sample size of sigma (if extrainfo=TRUE)}
#' \item{psi_eff_size}{Vector of effective sample sizes of psi (if extrainfo=TRUE)}
#' \item{xdist}{Vector of posterior draws from xdist between registered functions (if extrainfo=TRUE)}
#' \item{ydist}{Vector of posterior draws from ydist between registered functions (if extrainfo=TRUE)}
#' @references Lu, Y., Herbei, R., and Kurtek, S. (2017). Bayesian registration
#' of functions with a Gaussian process prior. Journal of Computational and
#' Graphical Statistics, DOI: 10.1080/10618600.2017.1336444.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @importFrom coda effectiveSize
#' @examples
#' data("simu_data")
#' myzpcn <- list(betas = c(0.1, 0.01, 0.005, 0.0001),
#'   probs = c(0.2, 0.2, 0.4, 0.2))
#' out = pair_align_functions_expomap(simu_data$f[,1], simu_data$f[,2],
#'   timet = simu_data$time, zpcn = myzpcn, extrainfo = TRUE)
#' # overall acceptance ratio
#' mean(out$accept)
#' # acceptance ratio by zpcn coefficient
#' with(out, tapply(accept, myzpcn$betas[betas.ind], mean))
pair_align_functions_expomap <- function(f1,
                                         f2,
                                         timet,
                                         iter = 2e4,
                                         burnin = min(5e3, iter / 2),
                                         alpha0 = 0.1,
                                         beta0 = 0.1,
                                         zpcn = list(betas = c(0.5, 0.05, 0.005, 0.0001),
                                                     probs = c(0.1, 0.1, 0.7, 0.1)),
                                         propvar = 1,
                                         init.coef = rep(0, 2 * 10),
                                         npoints = 200,
                                         extrainfo = FALSE) {
  # check input
  if (length(f1) != length(f2)) {
    stop('Length of f1 and f2 must be equal')
  }
  if (length(f1) != length(timet)) {
    stop('Length of f1 and timet must be equal')
  }
  if (length(zpcn$betas) != length(zpcn$probs)) {
    stop('In zpcn, betas must equal length of probs')
  }
  if ((length(init.coef) %% 2) != 0) {
    stop('Length of init.coef must be even')
  }

  # Number of sig figs to report in gamma_mat
  SIG_GAM <- 13

  # for now, back to list format of Yi's software
  f1 <- list(x = timet, y = f1)
  f2 <- list(x = timet, y = f2)

  # normalize timet to [0,1]
  # ([a,b] - a) / (b-a) = [0,1]
  rangex <- range(f1$x)
  f1$x <- (f1$x - rangex[1]) / (rangex[2] - rangex[1])
  f2$x <- (f2$x - rangex[1]) / (rangex[2] - rangex[1])

  # parameter settings
  pw.sim.global.burnin <- burnin
  valid.index <- pw.sim.global.burnin:iter
  pw.sim.global.Mg <- length(init.coef) / 2
  g.coef.ini <- init.coef
  numSimPoints <- npoints
  pw.sim.global.domain.par <-
    seq(from = 0,
        to = 1,
        length.out = numSimPoints)
  g.basis <- basis.fourier(f.domain = pw.sim.global.domain.par,
                           numBasis = pw.sim.global.Mg,
                           fourier.p = 1)
  sigma1.ini <- 1
  pw.sim.global.sigma.g <- propvar # proposal variance
  propose.g.coef <- function(g.coef.curr) {
    pCN.beta <- zpcn$betas
    pCN.prob <- zpcn$probs
    probm <- c(0, cumsum(pCN.prob))
    z = runif(1)
    for (i in 1:length(pCN.beta)) {
      if (z <= probm[i + 1] && z > probm[i]) {
        g.coef.new = rnorm(pw.sim.global.Mg * 2,
                           sd = pw.sim.global.sigma.g / rep(c(1:pw.sim.global.Mg), each =
                                                              2))
        result <-
          list(
            prop = sqrt(1 - pCN.beta[i] ^ 2) * g.coef.curr + pCN.beta[i] * g.coef.new,
            ind = i
          )
      }
    }
    return (result)
  }

  # srsf transformation
  q1 <- f.Q(f1)
  q2 <- f.Q(f2)

  # initialize chain
  obs.domain = q1$x

  if (min(f.exp1(
    f.basistofunction(
      f.domain = g.basis$x,
      coef = g.coef.ini,
      basis = g.basis
    )
  )$y) < 0) {
    stop("Invalid initial value of g")
  }

  ## results objects
  result.g.coef <- matrix(NA, ncol = iter, nrow = length(g.coef.ini))
  result.sigma1 <- rep(NA, iter)
  result.logl <- rep(NA, iter)
  result.SSE <- rep(NA, iter)
  result.accept <- rep(NA, iter)
  result.accept.betas <- rep(NA, iter)
  # matrix(0, nrow = length(zpcn$betas), ncol = 2,
  #  dimnames = list(beta = zpcn$betas, count = c('accept','attempt')))


  ## initialization
  g.coef.curr = g.coef.ini
  sigma1.curr = sigma1.ini
  SSE.curr = f.SSEg.pw(
    g = f.basistofunction(
      f.domain = g.basis$x,
      coef = g.coef.ini,
      basis = g.basis
    ),
    q1 = q1,
    q2 = q2
  )
  logl.curr = f.logl.pw(
    g = f.basistofunction(
      f.domain = g.basis$x,
      coef = g.coef.ini,
      basis = g.basis
    ),
    var1 = sigma1.ini ^ 2,
    q1 = q1,
    q2 = q2,
    SSEg = SSE.curr
  )

  result.g.coef[, 1] = g.coef.ini
  result.sigma1[1] = sigma1.ini
  result.SSE[1] = SSE.curr
  result.logl[1] = logl.curr

  ## update the chain for iter-1 times
  for (m in 2:iter) {
    #update g
    a = f.updateg.pw(
      g.coef.curr = g.coef.curr,
      g.basis = g.basis,
      var1.curr = sigma1.curr ^ 2,
      q1 = q1,
      q2 = q2,
      SSE.curr = SSE.curr,
      propose.g.coef = propose.g.coef
    )
    g.coef.curr = a$g.coef
    SSE.curr = a$SSE
    logl.curr = a$logl

    #update sigma1
    newshape = length(q1$x) / 2 + alpha0
    newscale = 1 / 2 * SSE.curr + beta0
    sigma1.curr = sqrt(rinvgamma(1, shape = newshape, scale = newscale))
    logl.curr = f.logl.pw(
      g = f.basistofunction(
        f.domain = g.basis$x,
        coef = g.coef.curr,
        basis = g.basis
      ),
      var1 = sigma1.curr ^ 2,
      q1 = q1,
      q2 = q2,
      SSEg = SSE.curr
    )

    #save update to results
    result.g.coef[, m] = g.coef.curr
    result.sigma1[m] = sigma1.curr
    result.SSE[m] = SSE.curr
    if (extrainfo) {
      result.logl[m] <- logl.curr
      result.accept[m] <- a$accept
      result.accept.betas[m] <- a$zpcnInd
    }
  }

  # calculate posterior mean of psi
  pw.sim.est.psi.matrix <-
    matrix(NA,
           nrow = length(pw.sim.global.domain.par),
           ncol = length(valid.index))
  for (k in 1:length(valid.index)) {
    g.temp = f.basistofunction(f.domain = pw.sim.global.domain.par,
                               coef = result.g.coef[, valid.index[k]],
                               basis = g.basis)
    psi.temp = f.exp1(g = g.temp)
    pw.sim.est.psi.matrix[, k] = psi.temp$y
  }

  result.posterior.psi.simDomain <- f.psimean(x = pw.sim.global.domain.par,
                                              y = pw.sim.est.psi.matrix)

  # resample to same number of points as the input f1 and f2.
  result.posterior.psi <- with(result.posterior.psi.simDomain,
                               approx(
                                 x = x,
                                 y = y,
                                 xout = f1$x,
                                 rule = 2
                               ))

  # transform posterior mean of psi to gamma
  ### aim to replace this by gamma mean below
  result.posterior.gamma <- f.phiinv(result.posterior.psi)

  # warped f2
  f2.warped <-
    with(result.posterior.gamma,
         warp_f_gamma(
           f = f2$y,
           time = x,
           gamma = y
         ))

  if (extrainfo) {
    # matrix of posterior draws from gamma
    gamma_mat <- apply(pw.sim.est.psi.matrix,
                       2,
                       function(vec) {
                         resamp <- approx(
                           x = pw.sim.global.domain.par,
                           y = vec,
                           xout = f1$x,
                           rule = 2
                         )
                         rawgam <- with(resamp, f.phiinv(list(x = x, y = y)))$y
                         round(rawgam, SIG_GAM)
                       })

    # x-distance, adapted from fdasrvf::elastic.distance
    ones <- rep(1, nrow(pw.sim.est.psi.matrix))
    Dx <- apply(pw.sim.est.psi.matrix,
                2,
                function(psi) {
                  v <- inv_exp_map(ones, psi)
                  sqrt(trapz(pw.sim.global.domain.par, v^2))
                })

    # y-distance, adapted from fdasrvf::elastic.distance
    Dy <- apply(gamma_mat,
                2,
                function(gam) {
                  q2warp <- warp_q_gamma(q2$y, q2$x, gam)
                  sqrt(trapz(q2$x, (q1$y - q2warp)^2))
                })

    # gamma stats: can add other stats of interest here...
    statsFun <- function(vec) {
      return(c(quantile(vec, probs = 0.025),
               quantile(vec, probs = 0.975)))
    }
    gamma.stats <- apply(gamma_mat, 1, statsFun)
    gamma_q025 <- gamma.stats['2.5%', ]
    gamma_q975 <- gamma.stats['97.5%', ]
  }

  # return object
  retVal <- list(
    f2_warped = f2.warped,
    gamma = result.posterior.gamma,
    g.coef = result.g.coef,
    psi = result.posterior.psi,
    sigma1 = result.sigma1
  )

  # optional return values
  if (extrainfo) {
    retVal$accept <- result.accept[-1]
    retVal$betas.ind <- result.accept.betas[-1]
    retVal$logl <- result.logl
    #    retVal$psi.mat <- pw.sim.est.psi.matrix
    retVal$gamma_mat <- gamma_mat
    retVal$gamma_q025 <- gamma_q025
    retVal$gamma_q975 <- gamma_q975
    # effective sample size
    sig.eff <- effectiveSize(result.sigma1)
    retVal$sigma_eff_size <- sig.eff
    psi.eff <- apply(pw.sim.est.psi.matrix, 1, effectiveSize)
    retVal$psi_eff_size <- psi.eff
    retVal$xdist <- Dx
    retVal$ydist <- Dy
  }

  return(retVal)
}
