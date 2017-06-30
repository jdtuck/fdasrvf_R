#' @export
summary.fdawarp <- function(object, ...){
  cat("-------------------------------------\n")
  cat(sprintf("Number of Functions: %d\n", ncol(object$f0)))
  cat(sprintf("Number of time points: %d\n", nrow(object$f0)))
  cat("-------------------------------------\n")
  cat(sprintf("Alignment with lambda: %f\n", object$lambda))
  cat(sprintf("Method: %s\n", object$method))
  cat(sprintf("Optimization Method: %s\n", object$omethod))
  cat(sprintf("Original Variance: %f\n", object$orig.var))
  cat(sprintf("Amplitude Variance: %f\n", object$amp.var))
  cat(sprintf("Phase Variance: %f\n", object$phase.var))
  cat(sprintf("Number of Iterations: %d\n", length(object$qun)-1))
}
