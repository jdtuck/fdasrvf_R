amplitude_distance <- function(q1, q2, scale = FALSE, norm_ratio = 1) {
  if (scale) {
    q1dotq2 <- innerprod_q2(q1, q2)
    if (q1dotq2 >  1) q1dotq2 <-  1
    if (q1dotq2 < -1) q1dotq2 <- -1
    d <- sqrt(acos(q1dotq2) ^ 2 + log(norm_ratio) ^ 2)
  } else {
    v <- q1 - q2
    d <- sqrt(innerprod_q2(v, v))
  }
  d
}

phase_distance <- function(gam) {
  M <- length(gam)
  grd <- seq(0, 1, length.out = M)
  binsize <- mean(diff(grd))
  psi <- sqrt(gradient(gam, binsize))
  v <- inv_exp_map(rep(1, M), psi)
  sqrt(trapz(grd, v ^ 2))
}

pointwise_mean <- function(q, beta, qmean,
                           basis = NULL,
                           scale = FALSE,
                           ms = "mean") {
  dims <- dim(q)
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]

  if (!scale) {
    # SVRF space has L2 geometry
    qmean <- rowMeans(q, dims = 2)
    betamean <- q_to_curve(qmean) + rowMeans(beta[, 1, ])
  } else {
    # SVRF space has L2 hypersphere geometry
    v <- v_d <- array(dim = c(L, M, N))
    d_i <- numeric(N)
    for (n in 1:N) {
      work_q <- q[, , n]
      q1dotq2 <- innerprod_q2(qmean, work_q)
      if (q1dotq2 >  1) q1dotq2 <-  1
      if (q1dotq2 < -1) q1dotq2 <- -1
      u <- work_q - q1dotq2 * qmean
      normu <- sqrt(innerprod_q2(u, u))
      if (normu > 1e-4)
        w <- u * acos(q1dotq2) / normu
      else
        w <- matrix(0, L, M)

      if (is.null(basis))
        v[, , n] <- w
      else
        v[, , n] <- project_tangent(w, qmean, basis)

      d_i[n] <- 0
      if (ms == "median") {
        # run for median only, saves computation time if getting mean
        d_i[n] <- sqrt(innerprod_q2(v[, , n], v[, , n])) # calculate sqrt of norm of v_i
        if (d_i[n] > 0)
          v_d[, , n] <- v[, , n] / d_i[n] # normalize v_i
        else
          v_d[, , n] <- v[, , n]
      } else
        v_d[, , n] <- matrix(0, nrow = L, ncol = M)
    }

    if (ms == "median") {
      # run for median only
      sumv <- rowSums(v_d, dims = 2)
      sum_dinv <- sum(1 / d_i)
      vbar <- sumv / sum_dinv
    } else {
      # run for mean only
      sumv <- rowSums(v, dims = 2)
      vbar <- sumv / N
    }

    normvbar <- sqrt(innerprod_q2(vbar, vbar))
    delta <- 0.5
    qmean <- cos(delta * normvbar) * qmean + sin(delta * normvbar) * vbar / normvbar
    if (!is.null(basis)) qmean <- project_curve(qmean)
    betamean <- q_to_curve(qmean)
  }

  list(qmean = qmean, betamean = betamean)
}
