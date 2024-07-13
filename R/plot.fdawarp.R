#' @importFrom graphics plot
#' @export
plot.fdawarp <- function(x, ...){

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

  M <- nrow(x$fn)
  mean_f0 <- rowMeans(x$f0)
  std_f0 <- apply(x$f0, 1, stats::sd)
  mean_fn <- rowMeans(x$fn)
  std_fn <- apply(x$fn, 1, stats::sd)

  graphics::matplot(
    x = (0:(M - 1)) / (M - 1),
    y = x$warping_functions,
    type = "l",
    col = colors,
    main = "Warping functions",
    xlab = "Time",
    ylab = "Warped Time"
  )

  graphics::matplot(
    x = x$time,
    y = x$fn,
    type = "l",
    col = colors,
    main = bquote(paste("Warped Data (", lambda == .(x$call$lambda), ")")),
    xlab ="Time",
    ylab = "Amplitude"
  )

  graphics::matplot(
    x = x$time,
    y = cbind(mean_f0, mean_f0 + std_f0, mean_f0 - std_f0),
    type = "l",
    lty = 1,
    col = c("black", "#66C2A5", "#66C2A5"),
    main = bquote(paste("Original Data: ", Mean %+-% STD)),
    xlab ="Time",
    ylab = "Amplitude"
  )

  graphics::matplot(
    x = x$time,
    y = cbind(mean_fn, mean_fn + std_fn, mean_fn - std_fn),
    type = "l",
    lty = 1,
    col = c("black", "#66C2A5", "#66C2A5"),
    main = bquote(paste(
      "Warped Data: ", Mean %+-% STD,
      " (", lambda == .(x$call$lambda), ")"
    )),
    xlab ="Time",
    ylab = "Amplitude"
  )

  if (x$call$centroid_type == "mean")
    plot(
      x = x$time,
      y = x$fmean,
      type = "l",
      col = "#66C2A5",
      main = bquote(paste(f[mean], " (", lambda == .(x$call$lambda), ")")),
      xlab ="Time"
      ,
      ylab = "Amplitude"
    )
  else
    plot(
      x = x$time,
      y = x$fmean,
      type = "l",
      col = "#66C2A5",
      main = bquote(paste(f[median], " (", lambda == .(x$call$lambda), ")")),
      xlab = "Time",
      ylab = "Amplitude"
    )
}
