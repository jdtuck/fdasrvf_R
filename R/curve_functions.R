
calculatecentroid <- function(beta,returnlength = F){
  n = nrow(beta)
  T1 = ncol(beta)

  betadot = apply(beta,1,gradient,1.0/(T1-1))
  betadot = t(betadot)

  normbetadot = apply(betadot,2,pvecnorm,2)
  integrand = matrix(0, n, T1)
  for (i in 1:T1){
    integrand[,i] = beta[,i] * normbetadot[i]
  }

  scale = trapz(seq(0,1,length.out=T1), normbetadot)
  centroid = apply(integrand,1,trapz,x = seq(0,1,length.out=T1))/scale
  if(returnlength)  return(list("length" = scale,"centroid" = centroid))
  return(centroid)
}

innerprod_q2 <- function(q1, q2) {
  T1 <- ncol(q1)
  sum(q1 * q2) / T1
}

find_best_rotation <- function(q1, q2) {
    eps <- .Machine$double.eps
    L <- nrow(q1)
    A <- q1 %*% t(q2)
    out <- svd(A)
    s <- out$d
    U <- out$u
    V <- out$v
    if (det(A) > 0)
      S <- diag(1, L)
    else {
      S <- diag(1, L)
      S[, L] <- -S[, L]
    }
    R <- U %*% S %*% t(V)
    q2new <- R %*% q2
    list(q2new = q2new, R = R)
}

calculate_variance <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    betadot = apply(beta,1,gradient,1.0/(T1-1))
    betadot = t(betadot)

    normbetadot = rep(0,T)
    centroid = calculatecentroid(beta)
    integrand = array(0, c(n,n,T1))
    time = seq(0,1,length.out=T1)
    for (i in 1:T1){
        normbetadot[i] = pvecnorm(betadot[,i],2)
        a1 = beta[,i] - centroid
        integrand[,,i] = a1 %*% t(a1) * normbetadot[i]
    }
    l = trapz(time, normbetadot)
    variance = trapz(time, integrand, 3)
    varaince = variance / l

    return(variance)
}


psi <- function(x, a, q){
    T1 = ncol(q)
    dim(a) = c(length(a),1)
    covmat = calculate_variance(x+repmat(a,1,T1))
    psi1 = covmat[1,1] - covmat[2,2]
    psi2 = covmat[1,2]
    psi3 = x[1,T1]
    psi4 = x[2,T1]

    return(list(psi1=psi1,psi2=psi2,psi3=psi3,psi4=psi4))
}


find_basis_normal <- function(q){
    n = nrow(q)
    T1 = ncol(q)

    f1 = matrix(0,n,T1)
    f2 = matrix(0,n,T1)
    for (i in 1:T1){
        f1[,i] = q[1,i]*q[,i]/pvecnorm(q[,i])+c(pvecnorm(q[,i]),0)
        f2[,i] = q[2,i]*q[,i]/pvecnorm(q[,i])+c(0,pvecnorm(q[,i]))
    }
    h3 = f1
    h4 = f2
    integrandb3 = rep(0,T1)
    integrandb4 = rep(0,T1)
    for (i in 1:T1){
        integrandb3[i] = t(q[,i])%*%h3[,i]
        integrandb4[i] = t(q[,i])%*%h4[,i]
    }
    b3 = h3 - q*trapz(seq(0,1,length.out=T1),integrandb3)
    b4 = h4 - q*trapz(seq(0,1,length.out=T1),integrandb4)

    basis = list(b3,b4)

    return(basis)
}


calc_j <- function(basis){
    b1 = basis[[1]]
    b2 = basis[[2]]
    T1 = ncol(b1)

    integrand11 = rep(0,T1)
    integrand12 = rep(0,T1)
    integrand22 = rep(0,T1)

    for (i in 1:T1){
        integrand11[i] = t(b1[,i])%*%b1[,i]
        integrand12[i] = t(b1[,i])%*%b2[,i]
        integrand12[i] = t(b2[,i])%*%b2[,i]
    }

    j = matrix(0,2,2)
    j[1,1] = trapz(seq(0,1,length.out=T1), integrand11)
    j[1,2] = trapz(seq(0,1,length.out=T1), integrand12)
    j[2,2] = trapz(seq(0,1,length.out=T1), integrand22)
    j[2,1] = j[1,2]

    return(j)
}


shift_f <- function(f, tau){
    n = nrow(f)
    T1 = ncol(f)
    fn = matrix(0, n, T1)
    if (tau == T1){
        fn[,(T1-tau+1):T1] = f[,1:tau]
    } else if (tau > 0){
        fn[,1:(T1-tau)] = f[,(tau+1):T1]
        fn[,(T1-tau+1):T1] = f[,1:tau]
    } else if (tau == 0) {
      fn = f
    } else {
        t = abs(tau)+1
        fn[,1:(T1-t+1)] = f[,(t):T1]
        fn[,(T1-t+2):T1] = f[,1:(t-1)]
    }
    return(fn)
}

find_rotation_seed_coord <- function(beta1, beta2,
                                     mode = "O",
                                     rotation = TRUE,
                                     scale = TRUE,
                                     lambda = 0.0) {
  L <- nrow(beta1)
  M <- ncol(beta1)
  out <- curve_to_q(beta1, scale = scale)
  q1 <- out$q
  len1 <- out$len

  scl <- 4
  minE <- Inf
  if (mode == "C")
    end_idx <- floor(M / scl)
  else
    end_idx <- 0

  for (ctr in 0:end_idx) {
    if (mode == "C") {
      if (scl * ctr <= end_idx){
        beta2n <- shift_f(beta2, scl * ctr)
        beta2n[, M] = beta2n[, 1]
      }
      else
        break
    } else
      beta2n <- beta2

    if (rotation) {
      out <- find_best_rotation(beta1, beta2n)
      beta2n <- out$q2new
      Rout <- out$R
    } else
      Rout <- diag(nrow(beta2n))

    out <- curve_to_q(beta2n, scale = scale)
    q2n <- out$q
    len2 <- out$len

    if (norm(q1 - q2n, 'F') > 0.0001) {
      if (scale) {
        q1i <- q1
        q2ni <- q2n
      } else {
        q1i <- q1 / sqrt(innerprod_q2(q1, q1))
        q2ni <- q2n / sqrt(innerprod_q2(q2n, q2n))
      }
      dim(q1i) <- M * L
      dim(q2ni) <- M * L
      gam0 <- DPQ(q1i, q2ni, L, M, lambda, 1, 0)
      gamI <- invertGamma(gam0)
      gam <- (gamI - gamI[1]) / (gamI[length(gamI)] - gamI[1])
      beta2new <- group_action_by_gamma_coord(beta2n, gam)
      out <- curve_to_q(beta2new, scale = scale)
      q2new <- out$q
      len2 <- out$len
      if (mode == "C") {
        q2new <- project_curve(q2new)
        beta2new <- q_to_curve(q2new)
      }
    } else {
      q2new <- q2n
      beta2new <- beta2n
      gam <- seq(0, 1, length.out = M)
    }

    Ec <- amplitude_distance(q1, q2new, scale = scale)

    if (Ec < minE) {
      Rbest <- Rout
      beta2best <- beta2new
      q2best <- q2new
      gambest <- gam
      minE <- Ec
      tau <- scl * ctr
    }
  }

  ratio <- 1
  if (scale)
    ratio <- len1 / len2

  list(
    beta2best = beta2best * ratio,
    q2best = q2best,
    Rbest = Rbest,
    gambest = gambest,
    tau = tau,
    scale = ratio
  )
}

find_rotation_seed_unique <- function(q1, q2,
                                      mode = "O",
                                      alignment = TRUE,
                                      rotation = TRUE,
                                      scale = TRUE,
                                      omethod = "DP",
                                      norm_ratio = 1.0,
                                      lambda = 0.0) {
  L <- nrow(q1)
  M <- ncol(q1)

  method <- match.arg(omethod, choices = c("DP", "DPo"))

  # Variables for DPQ2 algorithm
  grd <- seq(0, 1, length.out = M)
  nbhd_dim <- 7L

  scl <- 4
  minE <- Inf

  if (mode == "C")
    end_idx <- floor(M / scl)
  else
    end_idx <- 0

  for (ctr in 0:end_idx) {
    if (mode == "C")
      q2n <- shift_f(q2, scl * ctr)
    else
      q2n <- q2

    if (rotation) {
      out <- find_best_rotation(q1, q2n)
      q2n <- out$q2new
      Rbest <- out$R
    } else
      Rbest <- diag(nrow(q2n))

    if (norm(q1 - q2n, 'F') > 0.0001 && alignment) {
      if (scale) {
        q1i <- q1
        q2ni <- q2n
      } else {
        q1i <- q1 / sqrt(innerprod_q2(q1, q1))
        q2ni <- q2n / sqrt(innerprod_q2(q2n, q2n))
      }
      dim(q1i) <- M * L
      dim(q2ni) <- M * L

      switch(
        method,
        DP = {
          ret <- DPQ2(q1i, grd, q2ni, grd, L, M, M, grd, grd, M, M, lambda, nbhd_dim)
          Gvec <- ret$G[1:ret$size]
          Tvec <- ret$T[1:ret$size]
          gamI <- stats::approx(Tvec, Gvec, xout = grd)$y
        },
        DPo = {
          gamI <- DPQ(q2ni, q1i, L, M, lambda, 1, 0)
        }
      )

      gam <- (gamI - gamI[1]) / (gamI[length(gamI)] - gamI[1])
      q2new <- group_action_by_gamma(q2n, gam, scale = scale)

      if (mode == "C")
        q2new <- project_curve(q2new)
    } else {
      q2new <- q2n
      gam <- seq(0, 1, length.out = M)
    }

    Ec <- amplitude_distance(q1, q2new, scale = scale, norm_ratio = norm_ratio)

    if (Ec < minE) {
      Rbest <- Rbest
      q2best <- q2new
      gambest <- gam
      minE <- Ec
      tau <- scl * ctr
    }
  }

  list(
    q2best = q2best,
    Rbest = Rbest,
    gambest = gambest,
    tau = tau,
    d = minE
  )
}

find_rotation_and_seed_q <- function(q1,q2){
    n = nrow(q1)
    T1 = ncol(q1)
    Ltwo = rep(0,T1)
    Rlist = array(0,c(n,n,T1))
    for (ctr in 1:T1){
        q2n = shift_f(q2,ctr)
        out = find_best_rotation(q1,q2n)
        Ltwo[ctr] = innerprod_q2(q1-out$q2new,q1-out$q2new)
        Rlist[,,ctr] = out$R
    }

    tau = which.min(Ltwo)
    O_hat = Rlist[,,tau]
    q2new = shift_f(q2,tau)
    q2new = O_hat %*% q2new

    return(list(q2new=q2new,O_hat=O_hat,tau=tau))
}

group_action_by_gamma <- function(q, gamma, scale = TRUE) {
  L <- nrow(q)
  M <- ncol(q)
  grd <- seq(0, 1, length.out = M)
  gammadot <- gradient(gamma, 1.0 / M)
  qn <- matrix(nrow = L, ncol = M)

  for (l in 1:L)
    qn[l, ] <- stats::spline(grd, q[l, ], xout = gamma)$y * sqrt(gammadot)

  if (scale)
    qn <- qn / sqrt(innerprod_q2(qn, qn))

  qn
}

group_action_by_gamma_coord <- function(f, gamma) {
  L <- nrow(f)
  M <- ncol(f)
  fn <- matrix(nrow = L, ncol = M)
  grd <- seq(0, 1, length.out = M)

  for (l in 1:L)
    fn[l, ] <- stats::spline(grd, f[l, ], xout = gamma)$y

  fn
}


project_curve <- function(q){
    T1 = ncol(q)
    n = nrow(q)
    if (n==2){
      dt = 0.35
    }

    if (n==3) {
      dt = 0.2
    }

    epsilon = 1e-6

    e = diag(1,n)
    iter = 1
    res = rep(1,n)
    J = matrix(0,n,n)
    s = seq(0,1,length.out=T1)
    qnorm = rep(0,T1)
    G = rep(0,n)
    C = rep(0,301)

    qnew = q
    qnew = qnew / sqrt(innerprod_q2(qnew,qnew))

    while (pvecnorm(res) > epsilon){
        if (iter > 300){
            break
        }

        # Compute Jacobian
        for (i in 1:n){
            for (j in 1:n){
              J[i,j]  = 3 * trapz(s, qnew[i,]*qnew[j,])
            }
        }
        J = J + diag(1,n)

        for (i in 1:T1){
            qnorm[i] = pvecnorm(qnew[,i])
        }

        # Compute the residue
        for (i in 1:n){
            G[i] = trapz(s,qnew[i,]*qnorm)
        }

        res = -1*G

        if (pvecnorm(res)<epsilon)
            break

        x = solve(J,res)
        C[iter] = pvecnorm(res)

        delG = find_basis_normal(qnew)
        tmp = 0
        for (i in 1:n){
            tmp = tmp + x[i]*delG[[i]]*dt
        }
        qnew = qnew + tmp

        iter = iter + 1
    }

    qnew = qnew / sqrt(innerprod_q2(qnew,qnew))

    return(qnew)
}


pre_proc_curve <- function(beta, T1=100){
    beta = resamplecurve(beta,T1)
    q = curve_to_q(beta)$q
    qnew = project_curve(q)
    x = q_to_curve(qnew)
    a = -1*calculatecentroid(x)
    dim(a) = c(length(a),1)
    betanew = x + repmat(a,1,T1)
    A = diag(1,2)

    return(list(betanew=betanew,qnew=qnew,A=A))
}


inverse_exp_coord <- function(beta1, beta2, mode="O", rotated=T){
    T1 = ncol(beta1)
    centroid1 = calculatecentroid(beta1)
    dim(centroid1) = c(length(centroid1),1)
    beta1 = beta1 - repmat(centroid1, 1, T1)
    centroid2 = calculatecentroid(beta2)
    dim(centroid2) = c(length(centroid2),1)
    beta2 = beta2 - repmat(centroid2, 1, T1)

    q1 = curve_to_q(beta1)$q

    if (mode=="C"){
      isclosed = TRUE
    } else {
      isclosed = FALSE
    }

    # Iteratively optimize over SO(n) x Gamma using old DP
    out = reparam_curve(beta1, beta2, rotated=rotated, isclosed=isclosed, mode=mode)
    if (mode=="C")
      beta2 = shift_f(beta2, out$tau)

    beta2 = out$R %*% beta2
    gamI = invertGamma(out$gam)
    beta2 = group_action_by_gamma_coord(beta2, gamI)
    if (rotated){
      out = find_rotation_seed_coord(beta1, beta2, mode)
      q2n = curve_to_q(out$beta2new)$q
    } else {
      q2n = curve_to_q(beta2)$q
    }


    if (mode=="C"){
      q2n = project_curve(q2n)
    }

    # Compute geodesic distance
    q1dotq2 = innerprod_q2(q1/sqrt(innerprod_q2(q1, q1)), q2n/sqrt(innerprod_q2(q2n, q2n)))
    if (q1dotq2>1){
      q1dotq2 = 1.
    }

    dist = acos(q1dotq2)

    u = q2n - q1dotq2 * q1
    normu = sqrt(innerprod_q2(u,u))

    if (normu > 1e-4){
        v = u*acos(q1dotq2)/normu
    } else {
        v = matrix(0, nrow(beta1), T1)
    }

    return(list(v=v,dist=dist,gam=gamI))
}



inverse_exp <- function(q1, q2, beta2){
    T1 = ncol(q1)
    centroid1 = calculatecentroid(beta2)
    dim(centroid1) = c(length(centroid1),1)
    beta2 = beta2 - repmat(centroid1, 1, T1)

    # Optimize over SO(n) x Gamma
    beta1 = q_to_curve(q1)
    out = reparam_curve(beta1, beta2)
    gamI = invertGamma(out$gam)
    if (mode=="C")
      beta2 = shift_f(beta2, out$tau)

    beta2 = out$R %*% beta2

    # Applying optimal re-parameterization to the second curve
    beta2 = group_action_by_gamma_coord(beta2, gamI)
    q2 = curve_to_q(beta2)$q

    # Optimize over SO(n)
    out = find_rotation_and_seed_q(q1, q2)
    q2 = out$q2new

    # compute geodesic distance
    q1dotq2 = innerprod_q2(q1, q2)
    if (q1dotq2>1){
        q1dotq2 = 1.
    }

    dist = acos(q1dotq2)

    u = q2 - q1dotq2 * q1
    normu = sqrt(innerprod_q2(u,u))

    if (normu > 1e-4){
        v = u*acos(q1dotq2)/normu
    } else {
        v = matrix(0, 2, T1)
    }

    return(v)
}


gram_schmidt <- function(basis){
    b1 = basis[[1]]
    b2 = basis[[2]]

    basis1 = b1 / sqrt(innerprod_q2(b1,b1))
    b2 = b2 - innerprod_q2(basis1,b2)*basis1
    basis2 = b2 / sqrt(innerprod_q2(b2,b2))

    basis_o = list(basis1, basis2)

    return(basis_o)
}


project_tangent <- function(w, q, basis){
    w = w - innerprod_q2(w,q)*q
    bo = gram_schmidt(basis)

    wproj = w - innerprod_q2(w, bo[[1]])*bo[[1]] - innerprod_q2(w,bo[[2]])*bo[[2]]

    return(wproj)
}


scale_curve <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    normbetadot = rep(0,T1)
    betadot = matrix(0, n, T1)
    for (i in 1:n){
        betadot[i,] = gradient(beta[i,], 1.0/T1)
    }
    for (i in 1:T1){
        normbetadot[i] = pvecnorm(betadot[,i])
    }
    scale = trapz(seq(0,1,length.out=T1), normbetadot)
    beta_scaled = beta / scale

    return(list(beta_scaled=beta_scaled, scale=scale))
}


parallel_translate <- function(w, q1, q2, basis, mode='O'){
    wtilde = w - 2*innerprod_q2(w,q2) / innerprod_q2(q1+q2,q1+q2) * (q1+q2)
    l = sqrt(innerprod_q2(wtilde, wtilde))

    if (mode == 'C'){
        wbar = project_tangent(wtilde, q2, basis)
        normwbar = sqrt(innerprod_q2(wbar, wbar))
        if (normwbar>1e-4){
            wbar = wbar * l / normwbar
        }
    } else {
        wbar = wtilde
    }

    return(wbar)
}


rot_mat <- function(theta){
    O = matrix(0,2,2)
    O[1,1] = cos(theta)
    O[1,2] = -1*sin(theta)
    O[2,1] = sin(theta)
    O[2,2] = cos(theta)

    return(O)
}

karcher_calc <- function(q1, mu, basis,
                         rotated = TRUE,
                         scale = TRUE,
                         mode = "O",
                         lambda = 0.0,
                         ms = "mean") {
    tmp <- dim(q1)
    n <- tmp[1]
    T1 <- tmp[2]

    out <- find_rotation_seed_unique(
      q1 = mu,
      q2 = q1,
      mode = mode,
      rotation = rotated,
      scale = TRUE,
      lambda = lambda
    )

    qn_t <- out$q2best / sqrt(innerprod_q2(out$q2best, out$q2best))

    q1dotq2 <- innerprod_q2(mu, qn_t)

    if (q1dotq2 >  1) q1dotq2 <-  1
    if (q1dotq2 < -1) q1dotq2 <- -1

    dist <- acos(q1dotq2)

    u <- qn_t - q1dotq2 * q1
    normu <- sqrt(innerprod_q2(u, u))
    if (normu > 1e-4)
      w <- u * acos(q1dotq2) / normu
    else
      w <- matrix(0, n, T1)

    if (mode == "O")
      v <- w
    else
      v <- project_tangent(w, q1, basis)

    d_i <- 0
    if (ms == "median") {
      # run for median only, saves computation time if getting mean
      d_i <- sqrt(innerprod_q2(v, v)) # calculate sqrt of norm of v_i
      if (d_i > 0)
        v_d <- v / d_i #normalize v_i
      else
        v_d <- v
    } else
      v_d <- matrix(0, n, T1)

    list(v = v, v_d = v_d, dist = dist, d_i = d_i)
}

curve_align_sub <- function(beta1, q1, mu, mode, rotated, scale, lambda){
  out = find_rotation_seed_unique(mu, q1, mode, rotated, TRUE, lambda)
  gam = out$gambest
  rotmat = out$Rbest

  if (scale){
    q1n = out$q2best
    beta1n = q_to_curve(q1n)
  } else {
    beta1 = rotmat%*%beta1
    beta1n = group_action_by_gamma_coord(beta1, gam)
    q1n = curve_to_q(beta1n, scale)$q
  }

  T1 = ncol(beta1)
  centroid1 = calculatecentroid(beta1n)
  dim(centroid1) = c(length(centroid1), 1)
  beta1n = beta1n - repmat(centroid1, 1, T1)

  return(list(qn=q1n,betan=beta1n,rotmat=rotmat,gam=gam))
}


elastic_shooting <- function(q1, v, mode="O"){
    d = sqrt(innerprod_q2(v,v))
    if (d < 0.00001){
        q2n = q1
    } else {
        q2n = cos(d)*q1 + (sin(d)/d)*v
        if (mode == "C"){
        q2n = project_curve(q2n)
      }
    }

    return(q2n)
}

curve_to_srvf <- function(beta, scale = TRUE) {
  centroid <- calculatecentroid(beta)
  beta <- beta - centroid
  out <- curve_to_q(beta, scale = scale)
  list(q = out$q, qnorm = out$lenq, centroid = centroid)
}

match_f2_to_f1 <- function(srvf1, srvf2, beta2,
                           mode = "O",
                           alignment = TRUE,
                           rotation = FALSE,
                           scale = FALSE,
                           omethod = "DP",
                           include_length = FALSE,
                           lambda = 0.0) {
  norm_ratio <- 1
  if (include_length)
    norm_ratio <- srvf1$qnorm / srvf2$qnorm
  out <- find_rotation_seed_unique(
    srvf1$q, srvf2$q,
    mode = mode,
    alignment = alignment,
    rotation = rotation,
    scale = scale,
    lambda = lambda,
    omethod = omethod,
    norm_ratio = norm_ratio
  )

  # Compute amplitude distance
  d <- out$d

  # Compute phase distance
  gam <- out$gambest
  dx <- phase_distance(gam)

  # Compute beta2n
  qscale <- 1
  if (scale)
    qscale <- srvf1$qnorm / srvf2$qnorm
  betascale <- qscale^2

  if (mode == "C" && scale) {
    beta2n <- q_to_curve(out$q2best, scale = qscale)
  } else {
    beta2n <- beta2 - srvf2$centroid
    beta2n <- out$Rbest %*% beta2n
    beta2n <- group_action_by_gamma_coord(beta2n, gam)
    beta2n <- beta2n * betascale
  }
  beta2n <- beta2n + srvf1$centroid

  list(
    d = d,
    dx = dx,
    q2n = out$q2best,
    beta2n = beta2n,
    betascale = betascale,
    Rbest = out$Rbest,
    gambest = gam
  )
}

#' map shooting vector to curve at mean
#'
#'
#' @param v Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the shooting vectors
#' @param mu vector describing the mean
#' @param mode A character string specifying whether the input curves should be
#'   considered open (`mode == "O"`) or closed (`mode == "C"`). Defaults to
#'   `"O"`.
#' @param scale original scale of curve
#' @return A numeric array of the same shape as the input array `v` storing the
#'   curves of `v`.
#'
#' @keywords srvf alignment
#' @export
v_to_curve<-function(v, mu, mode="O", scale=1){
  n = nrow(mu)
  T1 = ncol(mu)
  if (ndims(v) == 0){
    dim(v) = c(n,T1)
    q2n = elastic_shooting(mu, v, mode)
    p = q_to_curve(q2n, scale)
  } else {

    p = matrix(0,T1,n)
    for (i in 1:n){
      v1 = v[,i]
      dim(v1) = c(n,T1)
      q2n = elastic_shooting(mu, v1, mode)
      p[,i] = q_to_curve(q2n, scale)
    }
  }
  return(p)
}
