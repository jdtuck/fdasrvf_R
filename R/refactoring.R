colSums_ext <- function(x, na.rm = FALSE, dims = 1) {
  if (is.null(dim(x))) return(x)
  colSums(x, na.rm = na.rm, dims = dims)
}

find_zero <- function(f, lower = 0, upper = 1, tol = 1e-8, maxiter = 1000L, ...) {
  # C_zeroin2 <- utils::getFromNamespace("C_zeroin2", "stats")
  # f.lower <- f(lower, ...)
  # f.upper <- f(upper, ...)
  # ff <- function(x) f(x, ...)
  # val <- .External2(C_zeroin2, ff, as.double(lower), as.double(upper), as.double(f.lower), as.double(f.upper), as.double(tol), as.integer(maxiter))
  val <- stats::uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter, ...)$root
  return(val[1])
}

integrate <- function(f, upper = 1) {
  stats::integrate(
    f = f,
    lower = 0,
    upper = upper,
    subdivisions = 10000L,
    stop.on.error = FALSE
  )$value
}

#' Converts a curve from matrix to functional data object
#'
#' @param beta A numeric matrix of size \eqn{L \times M} specifying a curve on
#'  an \eqn{L}-dimensional space observed on an evenly spaced grid of \eqn{[0,
#'  1]} of length \eqn{M}.
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'   1]} as input and returns the values of the original curve at \eqn{s}.
#' @export
#'
#' @examples
#' discrete2curve(beta[, , 1, 1])
discrete2curve <- function(beta) {
  dims <- dim(beta)
  L <- dims[1]
  M <- dims[2]
  grd <- seq(0, 1, length = M)
  funs <- lapply(1:L, \(l) {
    stats::splinefun(grd, beta[l, ], method = "natural")
  })
  \(s, deriv = 0) {
    out <- matrix(nrow = L, ncol = length(s))
    for (l in 1:L) {
      out[l, ] <- funs[[l]](s, deriv = deriv)
    }
    out
  }
}

#' Converts a warping function from vector to functional data object
#'
#' @param gam A numeric vector of length \eqn{M} specifying a warping function
#'   on an evenly spaced grid of \eqn{[0, 1]} of length \eqn{M}.
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'   1]} as input and returns the values of the original warping function at
#'   \eqn{s}.
#' @export
#'
#' @examples
#' discrete2warping(toy_warp$gam[, 1])
discrete2warping <- function(gam) {
  M <- length(gam)
  # step <- 1 / (M - 1)
  # left_slope <- (gam[2] - gam[1]) * step
  # right_slope <- (gam[M] - gam[M - 1]) * step
  # # y(x) = a (x - x0) + y0
  # gam <- c(gam[1] - step * left_slope , gam, gam[M] + step * right_slope)
  # grd <- seq(-step, 1 + step, length = M + 2)
  grd <- seq(0, 1, length = M)
  stats::splinefun(grd, gam, method = "hyman")
}

inverse_warping <- function(gamfun) {
  s <- seq(0, 1, length = 10000L)
  stats::splinefun(gamfun(s), s, method = "hyman")
  # \(s, deriv = 0) {
  #   if (deriv > 1)
  #     cli::cli_abort("The argument {.arg deriv} must be 0 or 1.")
  #   # Get gamma inverse
  #   gaminverse <- sapply(s, \(.s) {
  #     find_zero(\(t) gamfun(t) - .s, lower = 0, upper = 1, tol = 1e-8)
  #   })
  #   if (deriv == 0)
  #     return(gaminverse)
  #   1 / gamfun(gaminverse, deriv = 1)
  # }
}

#' Computes the centroid of a curve
#'
#' @param betafun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns its values at \eqn{s}.
#'
#' @return A numeric vector of length \eqn{L} storing the centroid of the curve.
#' @export
#'
#' @examples
#' betafun <- discrete2curve(beta[, , 1, 1])
#' get_curve_centroid(betafun)
get_curve_centroid <- function(betafun) {
  denom_integrand <- \(s) {
    betaprimevals <- betafun(s, deriv = 1)
    sqrt(apply(betaprimevals, 2, \(.x) sum(.x^2)))
  }
  denom_value <- integrate(denom_integrand)
  L <- dim(betafun(0))[1]
  num_values <- sapply(1:L, \(l) {
    num_integrand <- \(s) {
      betavals <- betafun(s)
      betaprimevals <- betafun(s, deriv = 1)
      normbetaprimevals <- sqrt(apply(betaprimevals, 2, \(.x) sum(.x^2)))
      betavals[l, ] * normbetaprimevals
    }
    integrate(num_integrand)
  })
  num_values / denom_value
}

#' Converts a curve to its SRVF representation
#'
#' @param beta A numeric matrix of size \eqn{L \times M} specifying a curve on
#'   an \eqn{L}-dimensional space observed on an evenly spaced grid of \eqn{[0,
#'   1]} of length \eqn{M}.
#' @param is_derivative A boolean value specifying whether the input \eqn{beta}
#'   is the derivative of the original curve. Defaults to `FALSE`.
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'   1]} as input and returns the values of the SRVF of the original curve at
#'   \eqn{s}.
#' @export
#'
#' @examples
#' curve2srvf(beta[, , 1, 1])
curve2srvf <- function(beta, is_derivative = FALSE) {
  if (!is_derivative) {
    if (is.matrix(beta))
      betafun <- discrete2curve(beta)
    else if (rlang::is_function(beta))
      betafun <- beta
    else
      cli::cli_abort("The argument {.arg beta} must be a matrix or a function.")
  } else {
    if (!rlang::is_function(beta))
      cli::cli_abort("The argument {.arg betaprime} must be a function")
    betaprime <- beta
  }

  \(s) {
    if (!is_derivative)
      fprime <- betafun(s, deriv = 1)
    else
      fprime <- betaprime(s)
    sqrt_fprime_norm <- sqrt(apply(fprime, 2, \(.x) sqrt(sum(.x^2))))
    is_zero <- sqrt_fprime_norm < sqrt(.Machine$double.eps)
    fprime[, is_zero] <- 0
    sqrt_fprime_norm[is_zero] <- 1
    fprime / matrix(
      data = sqrt_fprime_norm,
      nrow = nrow(fprime),
      ncol = length(s),
      byrow = TRUE
    )
  }
}

#' Converts from SRVF to curve representation
#'
#' @param qfun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the SRVF of an underlying
#'   curve at \eqn{s}.
#' @param beta0 A numeric vector of length \eqn{L} specifying the initial value
#'   of the underlying curve at \eqn{s = 0}.
#'
#' @return A function that takes a numeric vector \eqn{t} of values in \eqn{[0,
#'   1]} as input and returns the values of the underlying curve at \eqn{t}.
#' @export
#'
#' @examples
#' srvf2curve(curve2srvf(beta[, , 1, 1]))
srvf2curve <- function(qfun, beta0 = NULL) {
  L <- nrow(qfun(0))
  integrand <- \(s) {
    v <- qfun(s)
    v_norm <- sqrt(apply(v, 2, \(.x) sum(.x^2)))
    v * matrix(v_norm, nrow = L, ncol = length(s), byrow = TRUE)
  }
  Vectorize(\(t) {
    s <- seq(0, t, length = 10000L)
    integrand_matrix <- integrand(s)
    out <- trapz(s, integrand_matrix, dims = 2)
    if (!is.null(beta0))
      out <- out + beta0
    out
  })
}

#' Applies a warping function to a given curve
#'
#' @param betafun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the underlying curve at
#'   \eqn{s}.
#' @param gamfun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of a diffeomorphic warping
#'   function at \eqn{s}.
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'   1]} as input and returns the values of the warped curve at \eqn{s}.
#' @export
#'
#' @examples
#' curv <- discrete2curve(beta[, , 1, 1])
#' gamf <- discrete2warping(seq(0, 1, length = 100)^2)
#' warp_curve(curv, gamf)
warp_curve <- function(betafun, gamfun) {
  \(s) {
    betafun(gamfun(s))
  }
}

#' Applies a warping function to a given SRVF
#'
#' @param qfun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the SRVF of an underlying
#'   curve at \eqn{s}.
#' @param gamfun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of a diffeomorphic warping
#'   function at \eqn{s}.
#' @param betafun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the underlying curve at
#'   \eqn{s}. Defaults to `NULL`.
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'  1]} as input and returns the values of the warped SRVF.
#' @export
#'
#' @examples
#' q <- curve2srvf(beta[, , 1, 1])
#' warp_srvf(q, get_identity_warping())
warp_srvf <- function(qfun, gamfun, betafun = NULL) {
  L <- nrow(qfun(0))
  if (is.null(betafun)) {
    return(\(s) {
      gammaprime <- gamfun(s, deriv = 1)
      gammaprime[abs(gammaprime) < 1e-10] <- 0
      sqrt_gammaprime <- sqrt(gammaprime)
      sqrt_gammaprime <- matrix(
        data = sqrt_gammaprime,
        nrow = L,
        ncol = length(s),
        byrow = TRUE
      )
      qfun(gamfun(s)) * sqrt_gammaprime
    })
  }
  betaprime <- \(s) {
    gp_vals <- gamfun(s, deriv = 1)
    gp_vals <- matrix(gp_vals, nrow = L, ncol = length(s), byrow = TRUE)
    g_vals <- gamfun(s)
    betafun(g_vals, deriv = 1) * gp_vals
  }
  curve2srvf(betaprime, is_derivative = TRUE)
}

#' Computes the \eqn{L^2} distance between two SRVFs
#'
#' @param q1fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the first SRVF at \eqn{s}.
#' @param q2fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the second SRVF at \eqn{s}.
#' @param method A character string specifying the method to use for computing
#'   the \eqn{L^2} distance. Options are `"quadrature"` and `"trapz"`. Defaults
#'   to `"quadrature"`.
#'
#' @return A numeric value storing the \eqn{L^2} distance between the two SRVFs.
#' @export
#'
#' @examples
#' q1 <- curve2srvf(beta[, , 1, 1])
#' q2 <- curve2srvf(beta[, , 1, 2])
#' get_l2_distance(q1, q2)
get_l2_distance <- function(q1fun, q2fun, method = "quadrature") {
  if (method == "quadrature") {
    sqrt(integrate(\(s) colSums_ext((q1fun(s) - q2fun(s))^2)))
  } else if (method == "trapz") {
    grd <- seq(0, 1, length = 100000L)
    sqrt(trapz(grd, colSums_ext((q1fun(grd) - q2fun(grd))^2)))
  } else if (method == "simpson") {
    grd <- seq(0, 1, length = 100000L)
    sqrt(simpson(grd, colSums_ext((q1fun(grd) - q2fun(grd))^2)))
  } else if (method == "trapzCpp") {
    grd <- seq(0, 1, length = 100000L)
    sqrt(trapzCpp(grd, colSums_ext((q1fun(grd) - q2fun(grd))^2)))
  } else {
    cli::cli_abort("Invalid method")
  }
}

#' Computes the \eqn{L^2} inner product between two SRVFs
#'
#' @param q1fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the first SRVF at \eqn{s}.
#' @param q2fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the second SRVF at \eqn{s}.
#'
#' @return A numeric value storing the \eqn{L^2} inner product between the two
#'   SRVFs.
#' @export
#'
#' @examples
#' q1 <- curve2srvf(beta[, , 1, 1])
#' q2 <- curve2srvf(beta[, , 1, 2])
#' get_l2_inner_product(q1, q2)
get_l2_inner_product <- function(q1fun, q2fun) {
  integrate(\(s) colSums_ext(q1fun(s) * q2fun(s)))
}

#' Computes the \eqn{L^2} norm of an SRVF
#'
#' @param qfun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the SRVF at \eqn{s}.
#'
#' @return A numeric value storing the \eqn{L^2} norm of the SRVF.
#' @export
#'
#' @examples
#' q <- curve2srvf(beta[, , 1, 1])
#' get_l2_norm(q)
get_l2_norm <- function(qfun) {
  sqrt(get_l2_inner_product(qfun, qfun))
}

#' Projects an SRVF onto the Hilbert sphere
#'
#' @param qfun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the SRVF at \eqn{s}.
#' @param qnorm A numeric value specifying the \eqn{L^2} norm of the SRVF.
#'   Defaults to \eqn{\sqrt{\langle q, q \rangle}}.
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'   1]} as input and returns the values of the SRVF projected onto the Hilbert
#'   sphere.
#' @export
#'
#' @examples
#' q <- curve2srvf(beta[, , 1, 1])
#' to_hilbert_sphere(q)
to_hilbert_sphere <- function(qfun, qnorm = get_l2_norm(qfun)) {
  \(s) {
    q <- qfun(s)
    q / qnorm
  }
}

#' Computes the geodesic distance between two SRVFs on the Hilbert sphere
#'
#' @param q1fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the first SRVF at \eqn{s}.
#' @param q2fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the second SRVF at \eqn{s}.
#'
#' @return A numeric value storing the geodesic distance between the two SRVFs.
#' @export
#'
#' @examples
#' q1 <- curve2srvf(beta[, , 1, 1])
#' q2 <- curve2srvf(beta[, , 1, 2])
#' q1p <- to_hilbert_sphere(q1)
#' q2p <- to_hilbert_sphere(q2)
#' get_hilbert_sphere_distance(q1p, q2p)
get_hilbert_sphere_distance <- function(q1fun, q2fun) {
  ip_value <- get_l2_inner_product(q1fun, q2fun)
  if (ip_value >  1) ip_value <-  1
  if (ip_value < -1) ip_value <- -1
  acos(ip_value)
}

#' Computes the identity warping function
#'
#' @return A function that takes a numeric vector \eqn{s} of values in \eqn{[0,
#'   1]} as input and returns the values of the identity warping function at
#'   \eqn{s}.
#' @export
#'
#' @examples
#' get_identity_warping()
get_identity_warping <- function() {
  \(s, deriv = 0) {
    if (deriv > 1) return(rep(0, length(s)))
    if (deriv == 1) return(rep(1, length(s)))
    s
  }
}

#' Computes the distance between two warping functions
#'
#' @param gam1fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the first warping function
#'   at \eqn{s}.
#' @param gam2fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the second warping function
#'   at \eqn{s}.
#'
#' @return A numeric value storing the distance between the two warping
#'   functions.
#' @export
#'
#' @examples
#' gam1 <- discrete2warping(toy_warp$gam[, 1])
#' gam2 <- discrete2warping(toy_warp$gam[, 2])
#' get_warping_distance(gam1, gam2)
#' get_warping_distance(gam1, get_identity_warping())
get_warping_distance <- function(gam1fun, gam2fun) {
  psi1 <- \(s) {
    gam1prime <- gam1fun(s, deriv = 1)
    gam1prime[abs(gam1prime) < 1e-10] <- 0
    sqrt(gam1prime)
  }
  psi2 <- \(s) {
    gam2prime <- gam2fun(s, deriv = 1)
    gam2prime[abs(gam2prime) < 1e-10] <- 0
    sqrt(gam2prime)
  }
  theta <- get_hilbert_sphere_distance(psi1, psi2)
  if (theta < 1e-10) return(0)
  vfun <- \(s) theta / sin(theta) * (psi2(s) - cos(theta) * psi1(s))
  sqrt(get_l2_inner_product(vfun, vfun))
}

#' Computes the distance between two shapes
#'
#' @param q1fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the first SRVF at \eqn{s}.
#' @param q2fun A function that takes a numeric vector \eqn{s} of values in
#'   \eqn{[0, 1]} as input and returns the values of the second SRVF at \eqn{s}.
#' @param alignment A boolean value specifying whether the two SRVFs should be
#'   optimally aligned before computing the distance. Defaults to `FALSE`.
#' @param rotation A boolean value specifying whether the two SRVFs should be
#'   optimally rotated before computing the distance. Defaults to `FALSE`.
#' @param scale A boolean value specifying whether the two SRVFs should be
#'   projected onto the Hilbert sphere before computing the distance. Defaults
#'   to `FALSE`.
#' @param lambda A numeric value specifying the regularization parameter for the
#'   optimal alignment. Defaults to `0`.
#' @param M An integer value specifying the number of points to use when
#'   searching for the optimal rotation and alignment. Defaults to `200L`.
#'
#' @return A list with the following components:
#' - `amplitude_distance`: A numeric value storing the amplitude distance
#' between the two SRVFs.
#' - `phase_distance`: A numeric value storing the phase distance between the
#' two SRVFs.
#' - `optimal_warping`: A function that takes a numeric vector \eqn{s} of values
#' in \eqn{[0, 1]} as input and returns the optimal warping function evaluated
#' at \eqn{s}.
#' - `q1fun_modified`: A function that takes a numeric vector \eqn{s} of values
#' in \eqn{[0, 1]} as input and returns the first SRVF modified according to the
#' optimal alignment, rotation, and scaling.
#' - `q2fun_modified`: A function that takes a numeric vector \eqn{s} of values
#' in \eqn{[0, 1]} as input and returns the second SRVF modified according to
#' the optimal alignment, rotation, and scaling.
#' @export
#'
#' @examples
#' q1 <- curve2srvf(beta[, , 1, 1])
#' q2 <- curve2srvf(beta[, , 1, 2])
#' get_shape_distance(q1, q2)
get_shape_distance <- function(q1fun, q2fun,
                               alignment = FALSE,
                               rotation = FALSE,
                               scale = FALSE,
                               lambda = 0,
                               M = 200L) {
  L <- nrow(q1fun(0))
  if (nrow(q2fun(0)) != L)
    cli::cli_abort("The two input SRVFs should have the same dimension.")

  if (alignment || scale) {
    q1norm <- get_l2_norm(q1fun)
    q2norm <- get_l2_norm(q2fun)
  }

  if (scale) {
    q1fun_scaled <- to_hilbert_sphere(q1fun, q1norm)
    q2fun_scaled <- to_hilbert_sphere(q2fun, q2norm)
  } else {
    q1fun_scaled <- q1fun
    q2fun_scaled <- q2fun
  }

  if (rotation) {
    # Find optimal rotation
    grd <- seq(0, 1, length = M)
    R <- find_best_rotation(q1fun_scaled(grd), q2fun_scaled(grd))$R
    q2fun_scaled_rotated <- \(s) {R %*% q2fun_scaled(s)}
  } else {
    q2fun_scaled_rotated <- q2fun_scaled
  }

  if (alignment) {
    # Find optimal alignment
    grd <- seq(0, 1, length = M)

    if (scale) {
      Q1 <- q1fun_scaled(grd)
      Q2 <- q2fun_scaled_rotated(grd)
    } else {
      Q1 <- to_hilbert_sphere(q1fun_scaled, q1norm)(grd)
      Q2 <- to_hilbert_sphere(q2fun_scaled_rotated, q2norm)(grd)
    }
    dim(Q1) <- M * L
    dim(Q2) <- M * L

    # Variables for DPQ2 algorithm
    nbhd_dim <- 7L
    Gvec <- rep(0, M)
    Tvec <- rep(0, M)
    size <- 0

    ret <- DPQ2(Q1, grd, Q2, grd, L, M, M, grd, grd, M, M, Gvec, Tvec, size,
                lambda, nbhd_dim)

    Gvec <- ret$G[1:ret$size]
    Tvec <- ret$T[1:ret$size]
    # Numerical solution gamma_q2 to argmin d(q1, q2 o gamma)
    gamfun1 <- stats::splinefun(Tvec, Gvec, method = "hyman")
    # Numerical solution gamma_q1 to argmin d(q1 o gamma, q2)
    gamfun2 <- stats::splinefun(Gvec, Tvec, method = "hyman")
    # these two should be such that gamfun1(gamfun2^{-1}(s)) = s but it is not
    # the case because of numerical errors
    # Let's get the inverse of gamfun2 and we will pick between gamfun1 and
    # gamfun2_inv to minimize the distance and  get the optimal alignment of q2
    # with respect to q1.
    gamfun2_inv <- inverse_warping(gamfun2)

    work_q21 <- warp_srvf(q2fun_scaled_rotated, gamfun1) # q2 o gamfun1

    work_dist1 <- if (scale)
      get_hilbert_sphere_distance(q1fun_scaled, work_q21)
    else
      get_l2_distance(q1fun_scaled, work_q21)

    work_q22 <- warp_srvf(q2fun_scaled_rotated, gamfun2_inv) # q2 o gamfun2_inv

    work_dist2 <- if (scale)
      get_hilbert_sphere_distance(q1fun_scaled, work_q22)
    else
      get_l2_distance(q1fun_scaled, work_q22)

    if (work_dist1 < work_dist2) {
      gamfun <- gamfun1
      q2fun_scaled_rotated_aligned <- work_q21
      dist_amplitude <- work_dist1
    } else {
      gamfun <- gamfun2_inv
      q2fun_scaled_rotated_aligned <- work_q22
      dist_amplitude <- work_dist2
    }
  } else {
    gamfun <- get_identity_warping()
    q2fun_scaled_rotated_aligned <- q2fun_scaled_rotated
    dist_amplitude <- if (scale)
      get_hilbert_sphere_distance(q1fun_scaled, q2fun_scaled_rotated_aligned)
    else
      get_l2_distance(q1fun_scaled, q2fun_scaled_rotated_aligned)
  }

  dist_phase <- get_warping_distance(gamfun, get_identity_warping())

  list(
    ampitude_distance = dist_amplitude,
    phase_distance = dist_phase,
    optimal_warping = gamfun,
    q1fun_modified = q1fun_scaled,
    q2fun_modified = q2fun_scaled_rotated_aligned
  )
}

#' Computes the distance matrix between a set of shapes.
#'
#' @param qfuns A list of functions representing the shapes. Each function
#'   should be a SRVF. See [`curve2srvf()`] for more details on how to obtain
#'   the SRVF of a shape.
#' @inheritParams get_shape_distance
#' @param parallel_setup An integer value specifying the number of cores to use
#'   for parallel computing or an object of class `cluster`. In the latter case,
#'   it is used as is for parallel computing. If `parallel_setup` is `1L`, then
#'   no parallel computing will be used. The maximum number of cores to use is
#'   the number of available cores minus 1. Defaults to `1L`.
#'
#' @return An object of class `dist` containing the distance matrix between the
#'   shapes.
#' @export
#'
#' @examples
#' N <- 4L
#' srvfs <- lapply(1:N, \(n) curve2srvf(beta[, , 1, n]))
#' get_distance_matrix(srvfs, parallel_setup = 1L)
get_distance_matrix <- function(qfuns,
                                alignment = FALSE,
                                rotation = FALSE,
                                scale = FALSE,
                                lambda = 0,
                                M = 200L,
                                parallel_setup = 1L) {
  # Handles parallel computing strategy
  if (inherits(parallel_setup, "cluster")) {
    doParallel::registerDoParallel(parallel_setup)
    on.exit(parallel::stopCluster(parallel_setup))
  } else if (is.numeric(parallel_setup) && length(parallel_setup) == 1L) {
    ncores <- as.integer(parallel_setup)
    navail <- max(parallel::detectCores() - 1, 1)

    if (ncores > navail) {
      cli::cli_alert_warning(
        "The number of requested cores ({ncores}) is larger than the number of
      available cores ({navail}). Using the maximum number of available cores..."
      )
      ncores <- navail
    }

    if (ncores > 1L) {
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    } else
      foreach::registerDoSEQ()
  } else {
    cli::cli_abort("Invalid {.arg parallel_setup} argument. Must be a numeric
    scalar or an object of class {.cls cluster}.")
  }

  # Handles projection to Hilbert sphere if requested
  if (alignment || scale)
    qnorms <- sapply(qfuns, get_l2_norm)

  if (scale)
    qfuns_scaled <- mapply(to_hilbert_sphere, qfuns, qnorms, SIMPLIFY = FALSE)
  else
    qfuns_scaled <- qfuns

  L <- nrow(qfuns_scaled[[1]](0))
  N <- length(qfuns_scaled)
  K <- N * (N - 1) / 2

  k <- NULL
  out <- foreach::foreach(k = 0:(K - 1), .combine = cbind, .packages = "fdasrvf") %dopar% {
    # Compute indices i and j of distance matrix from linear index k
    i <- N - 2 - floor(sqrt(-8 * k + 4 * N * (N - 1) - 7) / 2.0 - 0.5)
    j <- k + i + 1 - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2

    # Increment indices as previous ones are 0-based while R expects 1-based
    q1fun_scaled <- qfuns[[i + 1]]
    q2fun_scaled <- qfuns[[j + 1]]

    if (rotation) {
      # Find optimal rotation
      grd <- seq(0, 1, length = M)
      R <- find_best_rotation(q1fun_scaled(grd), q2fun_scaled(grd))$R
      q2fun_scaled_rotated <- \(s) {R %*% q2fun_scaled(s)}
    } else {
      q2fun_scaled_rotated <- q2fun_scaled
    }

    if (alignment) {
      # Find optimal alignment
      grd <- seq(0, 1, length = M)

      if (scale) {
        Q1 <- q1fun_scaled(grd)
        Q2 <- q2fun_scaled_rotated(grd)
      } else {
        Q1 <- to_hilbert_sphere(q1fun_scaled, qnorms[i + 1])(grd)
        Q2 <- to_hilbert_sphere(q2fun_scaled_rotated, qnorms[j + 1])(grd)
      }
      dim(Q1) <- M * L
      dim(Q2) <- M * L

      # Variables for DPQ2 algorithm
      nbhd_dim <- 7L
      Gvec <- rep(0, M)
      Tvec <- rep(0, M)
      size <- 0

      ret <- DPQ2(Q1, grd, Q2, grd, L, M, M, grd, grd, M, M, Gvec, Tvec, size,
                  lambda, nbhd_dim)

      Gvec <- ret$G[1:ret$size]
      Tvec <- ret$T[1:ret$size]
      # Numerical solution gamma_q2 to argmin d(q1, q2 o gamma)
      gamfun1 <- stats::splinefun(Tvec, Gvec, method = "hyman")
      # Numerical solution gamma_q1 to argmin d(q1 o gamma, q2)
      gamfun2 <- stats::splinefun(Gvec, Tvec, method = "hyman")
      # these two should be such that gamfun1(gamfun2^{-1}(s)) = s but it is not
      # the case because of numerical errors
      # Let's get the inverse of gamfun2 and we will pick between gamfun1 and
      # gamfun2_inv to minimize the distance and  get the optimal alignment of q2
      # with respect to q1.
      gamfun2_inv <- inverse_warping(gamfun2)

      work_q21 <- warp_srvf(q2fun_scaled_rotated, gamfun1) # q2 o gamfun1

      work_dist1 <- if (scale)
        get_hilbert_sphere_distance(q1fun_scaled, work_q21)
      else
        get_l2_distance(q1fun_scaled, work_q21)

      work_q22 <- warp_srvf(q2fun_scaled_rotated, gamfun2_inv) # q2 o gamfun2_inv

      work_dist2 <- if (scale)
        get_hilbert_sphere_distance(q1fun_scaled, work_q22)
      else
        get_l2_distance(q1fun_scaled, work_q22)

      if (work_dist1 < work_dist2) {
        gamfun <- gamfun1
        q2fun_scaled_rotated_aligned <- work_q21
        dist_amplitude <- work_dist1
      } else {
        gamfun <- gamfun2_inv
        q2fun_scaled_rotated_aligned <- work_q22
        dist_amplitude <- work_dist2
      }
    } else {
      gamfun <- get_identity_warping()
      q2fun_scaled_rotated_aligned <- q2fun_scaled_rotated
      dist_amplitude <- if (scale)
        get_hilbert_sphere_distance(q1fun_scaled, q2fun_scaled_rotated_aligned)
      else
        get_l2_distance(q1fun_scaled, q2fun_scaled_rotated_aligned)
    }

    dist_phase <- get_warping_distance(gamfun, get_identity_warping())

    matrix(c(dist_amplitude, dist_phase), ncol = 1)
  }

  Da <- out[1, ]
  attributes(Da) <- NULL
  attr(Da, "Labels") <- 1:N
  attr(Da, "Size") <- N
  attr(Da, "Diag") <- FALSE
  attr(Da, "Upper") <- FALSE
  attr(Da, "call") <- match.call()
  attr(Da, "method") <- "amplitude"
  class(Da) <- "dist"

  Dp <- out[2, ]
  attributes(Dp) <- NULL
  attr(Dp, "Labels") <- 1:N
  attr(Dp, "Size") <- N
  attr(Dp, "Diag") <- FALSE
  attr(Dp, "Upper") <- FALSE
  attr(Dp, "call") <- match.call()
  attr(Dp, "method") <- "phase"
  class(Dp) <- "dist"

  list(Da = Da, Dp = Dp)
}
