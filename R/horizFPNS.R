#' Horizontal Functional Principal Component Analysis
#'
#' This function calculates vertical functional principal nested sphere analysis
#' on aligned data
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param var_exp compute no based on value percent variance explained (example: 0.95)
#' @param ci geodesic standard deviations (default = c(-1,0,1))
#' @param showplot show plots of principal directions (default = T)
#' @return Returns a hfpns object containing \item{gam_pns}{warping functions principal directions}
#' \item{psi_pns}{srvf principal directions}
#' \item{PNS}{PNS object}
#' @keywords srvf alignment
#' @references Yu, Q., Lu, X., and Marron, J. S. (2017), “Principal
#' Nested Spheres for Time-Warped Functional Data Analysis,” Journal
#' of Computational and Graphical Statistics, 26, 144–151.
#' @export
#' @examples
#' hfpns <- horizFPNS(simu_warp)
horizFPNS <- function(warp_data,
                      var_exp = 0.99,
                      ci = c(-1, 0, 1),
                      showplot = TRUE) {
  gam <- warp_data$warping_functions
  psi <- gam_to_psi(gam)
  TT <- nrow(psi)

  # PNS
  pnsdat <- apply(psi, 2, normalize_column)
  radius = mean(sqrt(apply(psi^2, 2, sum)))
  # get rough estimate of n.pc
  pca = eigen(stats::cov(t(pnsdat)))
  varExplained.psi = pca$values

  if (ncol(psi) < TT){
    n.pc = ncol(psi)
  } else {
    cs.psi = cumsum(varExplained.psi) / sum(varExplained.psi)
    n.pc = which(cs.psi >= 0.99)[1]
  }

  cli::cli_alert_info("Setting n.pc to {n.pc}...")

  obj = fastpns(pnsdat,
                n.pc = n.pc,
                sphere.type = "small",
                output = FALSE)

  varExplained = obj$percent
  cs = cumsum(obj$percent) / sum(obj$percent)

  no = which(cs >= var_exp)[1]

  proj_gam = array(0, dim = c(length(ci), nrow(psi), no))
  proj_psi = array(0, dim = c(length(ci), nrow(psi), no))
  for (j in 1:no) {
      std = stats::sd(obj$resmat[j, ])
      mean = mean(obj$resmat[j, ])
      dirtmp = ci * std + mean
      inmat = matrix(0, length(obj$PNS$radii), length(ci))
      inmat[j, ] = dirtmp
      PCvec = fastPNSe2s(inmat, obj)
      proj_psi[, , j] = as.matrix(PCvec) * radius
      for (k in 1:length(ci)){
        tmp = cumtrapz(seq(0, 1, length.out=TT), proj_psi[k, , j] * proj_psi[k, , j])
        proj_gam[k, , j] = (tmp - tmp[1]) / (tmp[length(tmp)] - tmp[1])
      }
  }

  hfpns = list()
  hfpns$gam_pns = proj_gam
  hfpns$psi_pns = proj_psi
  hfpns$no = no
  hfpns$psi = psi
  hfpns$PNS = obj$PNS
  hfpns$cumvar = cs
  hfpns$stds = ci
  hfpns$coef = obj$resmat
  hfpns$radius = radius
  hfpns$warp_data = warp_data

  class(hfpns) <- "hfpns"

  if (showplot) {
    plot(hfpns)
  }
  return(hfpns)
}
