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

pointwise_karcher_mean <- function(q, qmean,
                                   basis = NULL,
                                   scale = FALSE,
                                   ms = "mean",
                                   delta = 0.5) {
  dims <- dim(q)
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]

  v <- v_d <- array(dim = c(L, M, N))
  d_i <- numeric(N)
  for (n in 1:N) {
    w <- inverse_exponential_map(q[, , n], qmean, scale = scale)

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

  # Computes vbar
  if (ms == "median") {
    # run for median only
    sumv <- rowSums(v_d, dims = 2)
    sum_dinv <- sum(1 / d_i[d_i > 0])
    vbar <- sumv / sum_dinv
  } else {
    # run for mean only
    sumv <- rowSums(v, dims = 2)
    vbar <- sumv / N
  }

  qmean <- exponential_map(delta * vbar, qmean, scale = scale)

  if (!is.null(basis))
    qmean <- project_curve(qmean)

  list(qmean = qmean, vbar = vbar, v = v)
}

exponential_map <- function(v, f0, scale = FALSE) {
  if (scale) {
    normv <- sqrt(innerprod_q2(v, v))
    if (normv > sqrt(.Machine$double.eps))
      w <- v * sin(normv) / normv
    else
      w <- matrix(0, nrow = nrow(v), ncol = ncol(v))
    return(cos(normv) * f0 + w)
  }

  f0 + v
}

inverse_exponential_map <- function(f, f0, scale = FALSE) {
  if (scale) {
    cosval <- innerprod_q2(f0, f)
    if (cosval >  1) cosval <-  1
    if (cosval < -1) cosval <- -1
    theta <- acos(cosval)
    if (sin(theta) < sqrt(.Machine$double.eps))
      return(matrix(0, nrow = nrow(f), ncol = ncol(f)))
    else {
      v <- f - f0 * cosval
      return(theta * v / sin(theta))
    }
  }

  f - f0
}
