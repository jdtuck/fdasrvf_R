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
