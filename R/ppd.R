#' Compute Peak Persistence Diagram
#'
#' This computes the peak persistence diagram over a range of
#' lambda. This can help determine the proper elasticity (penalty).
#' This can be slow and recommended to run in parallel
#'
#' @param f A numeric matrix of shape \eqn{M \times N} specifying a sample of
#'   \eqn{N} curves observed on a grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the common grid on
#'   which all curves `f` have been observed.
#' @param max_lam maximum value of lambda. Defaults to `2`
#' @param num_lam number of steps. Defaults to `10`
#' @param pt the percentile of negative curvature of raw data Defaults to `.15`
#' @param penalty_method A string specifying the penalty term used in the
#'   formulation of the cost function to minimize for alignment. Choices are
#'   `"roughness"` which uses the norm of the second derivative, `"geodesic"`
#'   which uses the geodesic distance to the identity and `"norm"` which uses
#'   the Euclidean distance to the identity. Defaults to `"roughness"`.
#' @param centroid_type A string specifying the type of centroid to align to.
#'   Choices are `"mean"` or `"median"`. Defaults to `"mean"`.
#' @param center_warpings A boolean specifying whether to center the estimated
#'   warping functions. Defaults to `TRUE`.
#' @param smooth_data A boolean specifying whether to smooth curves using a box
#'   filter. Defaults to `FALSE`.
#' @param sparam An integer value specifying the number of times to apply the
#'   box filter. Defaults to `25L`. This is used only when `smooth_data = TRUE`.
#' @param parallel A boolean specifying whether to run calculations in parallel.
#'   Defaults to `TRUE`.
#' @param cores number of cores in parallel (default=-1, means all cores)
#' @param optim_method A string specifying the algorithm used for optimization.
#'   Choices are `"DP"`, `"DPo"`, and `"RBFGS"`. Defaults to `"DP"`.
#' @param max_iter An integer value specifying the maximum number of iterations.
#'   Defaults to `20L`.
#'
#' @return lam_opt optimal lam
#'
#' @keywords srsf alignment
#' @references Kim, Woo Min, Sutanoy Dasgupta, and Anuj Srivastava. "Peak-persistence diagrams for estimating shapes and functions from noisy data." arXiv preprint arXiv:2305.04826 (2023).
#'
#' @export
#' @examples
#' \dontrun{
#'   out <- ppd(simu_data$f, simu_data$time)
#' }
ppd <- function(f,
                time,
                max_lam = 2,
                num_lam = 10,
                pt = 0.15,
                penalty_method = c("roughness", "geodesic", "norm"),
                centroid_type = c("mean", "median"),
                center_warpings = TRUE,
                smooth_data = FALSE,
                sparam = 25L,
                parallel = TRUE,
                cores = -1,
                optim_method = c("DP", "DPo", "DP2", "RBFGS"),
                max_iter = 20L){

  penalty_method <- rlang::arg_match(penalty_method)
  centroid_type <- rlang::arg_match(centroid_type)
  optim_method <- rlang::arg_match(optim_method)

  # obtain a lambda candidate set
  lam_vec <- seq(0, max_lam, length.out = num_lam)
  fns <- list()
  for (i in 1:num_lam){
    obj <- time_warping(f, time, lam_vec[i],
                        penalty_method, centroid_type, center_warpings,
                        smooth_data, sparam, parallel, cores,
                        optim_method, max_iter)
    fns[[i]] <- obj$fn
  }

  # peak persistent diagram
  # get the trehsold for significant peak
  diff_t = mean(diff(time))
  taus = c()

  # compute tau values
  for (i in 1:ncol(f)){
    idx = findpeaks(f[, i])[,2]
    df2 = gradient(gradient(f[,i], diff_t), diff_t)
    tau = -df2 / max(-df2)
    tau[tau < 0] = 0
    taus = c(taus, c(tau[idx]))
  }
  th = stats::quantile(taus, pt)

  obj <- getPPDinfo(time, fns, lam_vec, th)
  persistent_peak_labels <- getPersistentPeaks(t(obj$IndicatorMatrix))

  # choose optimal lambda
  ref_row = rep(0, ncol(obj$IndicatorMatrix))
  ref_row[persistent_peak_labels] <- 1

  n_lams = nrow(obj$IndicatorMatrix)
  hamming_distances = rep(0, n_lams)
  exact_match_indices = c()

  for (i in 1:n_lams){
    comp = !is.nan(obj$IndicatorMatrix[i,])

    # check for exact match
    if (identical(comp, ref_row)){
      exact_match_indices = c(exact_match_indices, i)
    }
    # calculate Hamming distance
    hamming_distances[i] = sum(comp != ref_row)
  }

  if (length(exact_match_indices) > 0){
    idx_opt = min(exact_match_indices)
  } else {
    min_distance_indices = which.min(hamming_distances)
    idx_opt = min(min_distance_indices)
  }

  # draw PPD barchart
  drawPPDBarChart(obj$IndicatorMatrix, obj$Heights, lam_vec, idx_opt)
  # draw PPD surface
  drawPPDSurface(time, lam_vec, obj$FNm, obj$Heights, obj$Locs, obj$IndicatorMatrix, obj$Labels, idx_opt)

  return(lam_vec(idx_opt))

}

getPPDinfo <- function(t, Fa, lam, th){

  n_lams = length(lam)

  # compute the mean of each function in FN over its rows
  FNm = matrix(0, nrow=length(t), ncol=n_lams)
  for (i in 1:n_lams){
    FNm[,i] = rowMeans(Fa[[i]])
  }

  # find indices of lal maxima in the first function's mean
  idxMaxFirst = findpeaks(FNm[,1])[,2]

  # Initialize Labels and Locations for the first function
  Labels = list()
  Locs = list()
  Labels[[1]] = seq(1, length(idxMaxFirst))
  Locs[[1]] = idxMaxFirst

  # Initialize the maximum label number
  labelMax = max(Labels[[1]])

  # process each function to assign labels and locate peaks
  for (i in 1:(n_lams-1)){
    currentLabel = Labels[[i]]
    obj1 = peak_successor(Fa[[i]], Fa[[i+1]], currentLabel, labelMax)
    Labels[[i+1]] <- obj1$Labels
    labelMax = obj1$labelMax

    # find peak location is the next function's mean
    idxMaxNext = findpeaks(FNm[, i+1])[,2]
    Locs[[i+1]] <- idxMaxNext
  }

  # preprocess data to compute indicatorMatrix
  obj2 = PreprocessingForPPD(t, lam , Labels, Locs, labelMax, FNm, th)

  if (all(is.nan(obj2$Heights2))){
    warning('All peaks are ignored. A smaller threshold is required.')
  }

  obj = list()
  obj$IndicatorMatrix = obj2$IndicatorMatrix
  obj$Curvatures = obj2$Curvatures
  obj$Heights = obj2$Heights
  obj$Locs = Locs
  obj$Labels = Labels
  obj$FNm = FNm
  return(obj)

}

peak_successor <- function(f1, f2, labels1, labelMax){

  # combine f1 and f2 into a 3D array and compute the mean across the second dimension
  F <- array(c(f1, f2), dim = c(nrow(f1), ncol(f1), 2))
  fm <- apply(F, c(1, 3), mean)

  # compute peak ranges and labels for fm[,1]
  obj1 = computePeakRanges(fm[,1])
  ranges = obj1$ranges
  idx_max1 = obj1$idx_max

  # compute indices of local maxima in fm[,2]
  idx_max2 = findpeaks(fm[,2])[,2]

  if (length(idx_max1) == 0){
    labels2 = labelMax + seq(1,length(idx_max2))
  } else {
    # assign labels to peaks in fm[,2] based on matching ranges in fm[,1]
    labels2 = assignLabelsToPeaks(idx_max2, ranges, idx_max1, labels1)

    # ensure no overlapping labels in labels2
    labels2 = resolveOverlappingLabels(labels2, idx_max1, idx_max2, labels1)

    # assign new labels to unmatched peaks in fm[, 2]
    unmatched = labels2 == 0
    labels2[unmatched] = labelMax + seq(1,sum(unmatched))
    labelMax = labelMax + sum(unmatched)
  }

  obj = list()
  obj$Labels = labels2
  obj$labelMax = labelMax
  return(obj)
}

computePeakRanges <- function(data){
  # computes peak ranges defined by adjacent minima in the data
  idx_max = findpeaks(data)[,2]

  if (length(idx_max) == 0){
    warning("No peaks found in f1")
    obj = list()
    obj$ranges = c()
    obj$idx_max = c()
    return(obj)
  }

  idx_min = findpeaks(-1*data)[,2]
  idx_min = c(1, length(data), idx_min)
  idx_min = sort(unique(idx_min))

  ranges = matrix(0, nrow=length(idx_max), ncol=2)
  for (i in 1:length(idx_max)){
    ranges[i,1] = max(idx_min[idx_min<idx_max[i]])
    ranges[i,2] = min(idx_min[idx_min>idx_max[i]])
  }
  # remove degenerate ranges
  idx = which(ranges[,1] == ranges[,2])
  if (length(idx) > 0){
    ranges = ranges[-idx,]
  }

  obj = list()
  obj$ranges = ranges
  obj$idx_max = idx_max
  return(obj)
}

assignLabelsToPeaks <- function(idx_max2, ranges, idx_max1, labels1){
  # assigns labels to peaks in idx_max2 based on matching ranges in idx_max1
  labels = rep(0, length(idx_max2))
  for (i in 1:length(idx_max2)){
    # find the range in fm[,1] that contains the current peak in fm[,2]
    in_range = idx_max2[i] >= ranges[,1] & idx_max2[i] <= ranges[,2]
    matching_ranges = which(in_range)

    if (length(matching_ranges) > 0){
      # choose the cloest peak if multiple ranges match
      if (length(matching_ranges) > 1){
        closest_idx = which.min(abs(idx_max1[matching_ranges] - idx_max2[i]))
        matching_range = matching_ranges[closest_idx]
      } else {
        matching_range = matching_ranges
      }
      labels[i] = labels1[matching_range]
    }
  }

  return(labels)

}

resolveOverlappingLabels <- function(labels, idx_max1, idx_max2, labels1){
  # ensures no overlapping labels in the assigned labels
  unique_labels = unique(labels[labels > 0])
  for (label in unique_labels){
    duplicates = which(labels==label)
    if (length(duplicates) > 1){
      # keep the cloest peak and reset others
      distances = abs(idx_max1[labels==label] - idx_max2[duplicates])
      min_idx = which.min(distances)
      duplicates = duplicates[-min_idx]
      labels[duplicates] = 0
    }
  }
  return(labels)
}

PreprocessingForPPD <- function(t, lam, Labels, Locs, labelMax, FNm, th){
  K = length(lam)

  # Initialize output matrices
  IndicatorMatrix = matrix(NaN, nrow=K, ncol=labelMax)
  curvatures = matrix(0, nrow=K, ncol=labelMax)
  heights = matrix(NaN, nrow=K, ncol=labelMax)

  # assume t is uniformly spaced; compute time step
  dx = t[2] - t[1]

  # process each function
  for (i in 1:K){
    # extract the function values for the current parameter lambda
    fnm = FNm[, i]

    # compute negative curvature (second derivative)
    negCurvature = -gradient(gradient(fnm, dx), dx)

    # ensure non-negative curvature values
    negCurvature[negCurvature < 0] <- 0

    # normalize non-negative curvature to [0,1] if possible
    maxNegCurvature = max(negCurvature)
    if (maxNegCurvature > 0){
      negCurvature = negCurvature / maxNegCurvature
    }

    # retrieve peak locations and labels for the current function
    locsCurrent = Locs[[i]]
    labelsCurrent = Labels[[i]]

    # select negative curvature values at specified peak locations
    negCurvSelected = negCurvature[locsCurrent]

    # update curvatures and heights matrices at the appropriate labels
    curvatures[i, labelsCurrent] = negCurvSelected
    heights[i, labelsCurrent] = fnm[locsCurrent]

    # apply threshold to select significant peaks based on curvature
    significantLabels = labelsCurrent[negCurvSelected >= th]

    # update the indicator matrix for significant peaks
    IndicatorMatrix[i, significantLabels] =1
  }

  # compute heights 2 by multiplying heights with the indicator matrix
  heights2 = IndicatorMatrix * heights

  obj = list()
  obj$IndicatorMatrix = IndicatorMatrix
  obj$curvatures = curvatures
  obj$Heights = heights
  obj$Heights2 = heights2
  return(obj)
}

getPersistentPeaks <- function(IndicatorMatrix){
  if (length(IndicatorMatrix) == 0){
    Clt2 = c()
    return(Clt2)
  }

  # count the number of ones (occurrences) for each peak (ignore NaNs)
  occurrenceCounts = rowSums(IndicatorMatrix, na.rm = TRUE)

  data = c(occurrenceCounts, 0)
  if (all(data==0) | length(unique(occurrenceCounts)) == 1){
    Clt2 = c()
    return(Clt2)
  }

  # compute pairwise distances between observations
  pairwiseDistances <- stats::dist(data)
  hc <- stats::hclust(pairwiseDistances, method = "ward.D2")
  clusterAssignments = hc$labels2
  referenceCluster = clusterAssignments[length(clusterAssignments)]
  clusterAssignments = clusterAssignments[-length(clusterAssignments)]

  # identify indices where cluster assignments differ from the reference
  Clt2 = which(clusterAssignments != referenceCluster)
  return(Clt2)

}

drawPPDBarChart <- function(IndicatorMatrix, Heights, lam, idx_opt){
  lam_diff = lam[2] - lam[1]
  len_lam = nrow(IndicatorMatrix)
  labelMax = ncol(IndicatorMatrix)

  plot(c(lam[1], lam[length(lam)]), c(0.5, labelMax+0.5), type = "n", xlab = "lambda", ylab = "Peak Index",
       main = "", yaxt='n')

  for (i in 1:len_lam){
    label_all_peaks = which(!is.nan(Heights[i,]))
    label_persistent_peaks = which(!is.nan(IndicatorMatrix[i,]))

    if (sum(label_all_peaks) == 0){
      next
    }

    for (j in 1:length(label_all_peaks)){
      x = lam[i]
      y = label_all_peaks[j]

      x1 = x - lam_diff/2
      x2 = x + lam_diff/2
      y1 = y - 0.5
      y2 = y + 0.5

      if (y %in% label_persistent_peaks){
        graphics::rect(x1, y1, x2, y2, col="black")
      } else {
        graphics::rect(x1, y1, x2, y2, col="#717171")
      }
    }
  }

  for (j in 1:labelMax){
    graphics::abline(h = j+0.5, col = "blue", lwd = .5)
  }

  graphics::abline(v = lam[idx_opt], col = "magenta", lty=2, lwd = 2)

  ticks = seq(1,labelMax)
  graphics::axis(side = 2, at = ticks)

}

drawPPDSurface <- function(t,lam,FNm,Heights,Locs,IndicatorMatrix,Labels,idx_opt){
  n_lams = length(lam)
  labelMax = ncol(IndicatorMatrix)

  LocationMatrix_full = matrix(NaN, nrow=n_lams, ncol=labelMax)
  for (i in 1:n_lams){
    LocationMatrix_full[i, Labels[[i]]] = Locs[[i]]
  }
  LocationMatrix_sig = LocationMatrix_full * IndicatorMatrix
  HeightMatrix_full = Heights
  HeightMatrix_sig = HeightMatrix_full * IndicatorMatrix

  if (requireNamespace("plot3Drgl", quietly = TRUE)){
    plot3D::persp3D(
      x = t,
      y = lam,
      z = t(FNm),
      col = viridisLite::viridis(128),
      plot = FALSE,
      main = "",
      xlab="t",
      ylab="g_lambda",
      zlab="lambda",
      ticktype = "detailed",
      box = FALSE
    ) +
      plot3D::lines3D(
        x = t,
        y = lam[idx_opt]*rep(1,length(t)),
        z = FNm[, idx_opt],
        col = "magenta",
        lty=2,
        lwd = 3,
        add = TRUE,
        plot = FALSE
      )
    for (j in 1:labelMax){
      # find non-NaN indices for full location matrix and plot
      idx_full = which(!is.nan(LocationMatrix_full[,j]))
      plot3D::points3D(
        x = t[LocationMatrix_full[idx_full, j]],
        y = lam[idx_full],
        z = HeightMatrix_full[idx_full, j],
        col = "black",
        lty = 1,
        lwd = 1.5,
        add = TRUE,
        plot = FALSE
      )

      # find non-naN indices for significant location matrix and plot
      idx_sig = which(!is.nan(LocationMatrix_sig[, j]))
      plot3D::points3D(
        x = t[LocationMatrix_sig[,j]],
        y = lam[idx_sig],
        z = HeightMatrix_sig[idx_sig, j],
        col = "black",
        lty = 2,
        lwd = 2,
        add = TRUE,
        plot = FALSE
      )

    }

    plot3Drgl::plotrgl()
    rgl::par3d("windowRect" = c(0, 0, 640, 640))
    rgl::grid3d(c("x", "y+", "z"))
    rgl::axes3d(c('x--', "y--", 'z'))
    rgl::title3d(xlab = "t", ylab = "g_lambda", zlab="lamda")
  } else {
    graphics::image(t, lam, t(FNm), main = "", xlab="t", zlab="g_lambda",
                    ylab = "lambda", col = viridisLite::viridis(128))
    graphics::lines(t, lam(idx_opt)*rep(1, length(t)), col="magenta", lwd=2)

    for (j in 1:labelMax){
      # find non-NaN indices for full location matrix and plot
      idx_full = which(!is.nan(LocationMatrix_full[,j]))
      graphics::points(t[LocationMatrix_full[idx_full, j]], lam[idx_full],
                       col="black", lwd=1.5)

      # find non-naN indices for significant location matrix and plot
      idx_sig = which(!is.nan(LocationMatrix_sig[, j]))
      graphics::points(t[LocationMatrix_sig[,j]], lam[idx_sig], col="black",
                       lwd = 2)

    }
  }

}
