#' Plot functional data
#'
#' This function plots functional data on a time grid
#'
#' @param t A numeric vector of length \eqn{M} specifying the common grid on
#'   which all curves `f` have been observed.
#' @param f A numeric matrix of shape \eqn{M \times N} specifying a sample of
#'   \eqn{N} curves observed on a grid of size \eqn{M}.
#'
#' @importFrom graphics plot
#' @export
f_plot <- function(t, f, ...){

  colors = c(
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#B3B3B3"
  )

  M <- nrow(f)

  graphics::matplot(
    x = t,
    y = f,
    type = "l",
    col = colors,
    main = "Original Data",
    xlab = "Time",
    ylab = "Amplitude"
  )
}
