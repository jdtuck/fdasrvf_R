#' Elastic Functional Data Analysis
#'
#' A library for functional data analysis using the square root velocity
#' framework which performs pair-wise and group-wise alignment as well as
#' modeling using functional component analysis.
#'
#' @name fdasrvf
#'
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using Fisher-Rao metric,
#'   arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using phase and amplitude separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references J. D. Tucker, W. Wu, and A. Srivastava, Phase-amplitude
#'   separation of proteomics data using extended Fisher-Rao metric, Electronic
#'   Journal of Statistics, Vol 8, no. 2. pp 1724-1733, 2014.
#' @references J. D. Tucker, W. Wu, and A. Srivastava, ``Analysis of signals
#'   under compositional noise with applications to SONAR data," IEEE Journal of
#'   Oceanic Engineering, Vol 29, no. 2. pp 318-330, Apr 2014.
#' @references Tucker, J. D. 2014, Functional Component Analysis and Regression
#'   using Elastic Methods. Ph.D. Thesis, Florida State University.
#' @references Robinson, D. T. 2012, Function Data Analysis and Partial Shape
#'   Matching in the Square Root Velocity Framework. Ph.D. Thesis, Florida State
#'   University.
#' @references Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds
#'   with Applications. Ph.D. Thesis, Florida State University.
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian
#'   registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'   analysis of elastic curves in euclidean spaces. Pattern Analysis and
#'   Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian
#'   registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @references W. Xie, S. Kurtek, K. Bharath, and Y. Sun, A geometric approach
#'   to visualization of variability in functional data, Journal of American
#'   Statistical Association 112 (2017), pp. 979-993.
#' @references Lu, Y., R. Herbei, and S. Kurtek, 2017: Bayesian registration of
#'   functions with a Gaussian process prior. Journal of Computational and
#'   Graphical Statistics, 26, no. 4, 894–904.
#' @references Lee, S. and S. Jung, 2017: Combined analysis of amplitude and
#'   phase variations in functional data. arXiv:1603.01775, 1–21.
#' @references J. D. Tucker, J. R. Lewis, and A. Srivastava, “Elastic Functional
#'   Principal Component Regression,” Statistical Analysis and Data Mining, vol.
#'   12, no. 2, pp. 101-115, 2019.
#' @references J. D. Tucker, J. R. Lewis, C. King, and S. Kurtek, “A Geometric
#'   Approach for Computing Tolerance Bounds for Elastic Functional Data,”
#'   Journal of Applied Statistics, 10.1080/02664763.2019.1645818, 2019.
#' @references T. Harris, J. D. Tucker, B. Li, and L. Shand, "Elastic depths for
#'   detecting shape anomalies in functional data," Technometrics,
#'   10.1080/00401706.2020.1811156, 2020.
#' @references J. D. Tucker and D. Yarger, “Elastic Functional Changepoint
#'   Detection of Climate Impacts from Localized Sources”, Envirometrics,
#'   10.1002/env.2826, 2023.
#'
#' @docType package
#' @useDynLib fdasrvf, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %dopar%
#' @aliases fdasrvf fdasrvf-package
NULL
