#' @export
summary.fdawarp <- function(out){
  cat("-------------------------------------\n")
  cat(sprintf("Number of Functions: %d\n", ncol(out$f0)))
  cat(sprintf("Number of time points: %d\n", nrow(out$f0)))
  cat("-------------------------------------\n")
  cat(sprintf("Alignment with lambda: %f\n", out$lambda))
  cat(sprintf("Method: %s\n", out$method))
  cat(sprintf("Optimization Method: %s\n", out$omethod))
  cat(sprintf("Original Variance: %f\n", out$orig.var))
  cat(sprintf("Amplitude Variance: %f\n", out$amp.var))
  cat(sprintf("Phase Variance: %f\n", out$phase.var))
  cat(sprintf("Number of Iterations: %d\n", length(out$qun)-1))
}
