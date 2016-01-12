#' Elastic Shape Distance
#'
#' Calculate elastic shape distance between two curves beta1 and beta2
#'
#' @param beta1 array describing curve1 (n,T)
#' @param beta2 array describing curve
#' @return d geodesic distance
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("curve_data")
#' d = calc_shape_dist(curve_data$beta1,curve_data$beta2)
calc_shape_dist <- function(beta1, beta2){
    out = inverse_exp_coord(beta1, beta2)

    return(out$d)
}
