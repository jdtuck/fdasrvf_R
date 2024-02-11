#' Plot Curve
#'
#' This function plots open or closed curves
#'
#' @param beta Array of sizes \eqn{n \times T} describing a
#'     curve of dimension \eqn{n} evaluated on \eqn{T} points
#' @param add add to current plot (default = `TRUE`)
#' @param ... additional plotting parameters
#' @return Return shape confidence intervals
#' @keywords bootstrap
#' @export
plot_curve <- function(beta, add=FALSE, ...){

  if (add){
    lines(beta[1,], beta[2,], ...)
  } else {
    plot(beta[1,], beta[2,], type="l", ...)
  }
}
