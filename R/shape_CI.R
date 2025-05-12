#' Shape Confidence Interval Calculation using Bootstrap Sampling
#'
#' This function computes Confidence bounds for shapes using elastic metric
#'
#' @param beta Array of sizes \eqn{n \times T \times N} describing \eqn{N}
#'     curves of dimension \eqn{n} evaluated on \eqn{T} points
#' @param a confidence level (default = 0.95)
#' @param no number of principal components (default = 5)
#' @param Nsamp number of functions to generate (default = 100)
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotated Optimize over rotation (default = `TRUE`)
#' @param scale scale curves to unit length (default = `TRUE`)
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param parallel enable parallel processing (default = T)
#' @return Return shape confidence intervals
#' @keywords bootstrap
#' @concept bounds
#' @references J. D. Tucker, J. R. Lewis, C. King, and S. Kurtek, “A Geometric
#'   Approach for Computing Tolerance Bounds for Elastic Functional Data,”
#'   Journal of Applied Statistics, 10.1080/02664763.2019.1645818, 2019.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative Models for
#'   Function Data using Phase and Amplitude Separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
shape_CI <- function(beta, a=.95, no = 5, Nsamp=100,
                     mode = "O", rotated = TRUE, scale = TRUE,
                     lambda = 0.0, parallel=TRUE){

  # Align Data --------------------------------------------------------------
  out.med <- multivariate_karcher_mean(beta, mode=mode, rotation=rotated, scale=scale,
                                       lambda=lambda, ms="median")

  samples <- sample_shapes(out.med, no, Nsamp)
  box <- curve_boxplot(samples, alpha = 1 - a)

  return(list(Q1a=box$Q1a,Q3a=box$Q3a, samples=samples, out.med=out.med, box=box))

}
