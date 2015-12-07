#' Find optimum reparameterization between two images
#'
#' This function aligns to images using the q-map framework
#'
#' @param It template image matrix
#' @param Im test image matrix
#' @param gam initial warping array
#' @param b basis matrix
#' @param stepsize gradient stepsize (default=1e-5)
#' @param itermax maximum number of iterations (default=1000)
#' @param lmark use landmarks (default=FALSE)
#' @return Returns a list containing \item{gamnew}{final warping}
#' \item{Inew}{aligned image}
#' \item{H}{energy}
#' \item{stepsize}{final stepsize}
#' @keywords image alignment
#' @references Q. Xie, S. Kurtek, E. Klassen, G. E. Christensen and A. Srivastava. Metric-based pairwise and multiple image registration. IEEE European Conference on Computer Vision (ECCV), September, 2014
#' @export
reparam_image <- function(It, Im, gam, b, stepsize=1e-5, itermax=1000, lmark=FALSE){
    
}
