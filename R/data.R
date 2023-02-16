#' Simulated two Gaussian Dataset
#'
#' A functional dataset where the individual functions are given by: \eqn{y_i(t)
#' = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3],
#' ~i=1,2,\dots, 21}, where \eqn{z_{i,1}} and \eqn{z_{i,2}} are *i.i.d.* normal
#' with mean one and standard deviation 0.25. Each of these functions is then
#' warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} -
#' 1}) - 3} if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}}
#' (\eqn{gamma_{id}(t) = t}) is the identity warping). The variables are as
#' follows: f containing the 21 functions of 101 samples and time which
#' describes the sampling.
#'
#' @format ## `simu_data`
#' A list with 2 components:
#'
#' - `f`: A numeric matrix of shape \eqn{101 \times 21} storing a sample of size
#' \eqn{N = 21} of curves evaluated on a grid of size \eqn{M = 101}.
#' - `time`: A numeric vector of size \eqn{M = 101} storing the grid on which
#' the curves `f` have been evaluated.
#'
#' @source ???
"simu_data"

#' Distributed Gaussian Peak Dataset
#'
#' A functional dataset where the individual functions are given by a Gaussian
#' peak with locations along the \eqn{x}-axis. The variables are as follows: f
#' containing the 29 functions of 101 samples and time which describes the
#' sampling.
#'
#' @format ## `toy_data`
#' A list with two components:
#'
#' - `f`: A numeric matrix of shape \eqn{101 \times 29} storing a sample of size
#' \eqn{N = 29} of curves evaluated on a grid of size \eqn{M = 101}.
#' - `time`: A numeric vector of size \eqn{M = 101} storing the grid on which
#' the curves `f` have been evaluated.
#'
#' @source ???
"toy_data"

#' Berkeley Growth Velocity Dataset
#'
#' Combination of both boys and girls growth velocity from the Berkley dataset.
#'
#' @format ## `growth_vel`
#' A list with two components:
#'
#' - `f`: A numeric matrix of shape \eqn{69 \times 93} storing a sample of size
#' \eqn{N = 93} of curves evaluated on a grid of size \eqn{M = 69}.
#' - `time`: A numeric vector of size \eqn{M = 69} storing the grid on which
#' the curves `f` have been evaluated.
#'
#' @source ???
"growth_vel"

#' Aligned Distributed Gaussian Peak Dataset
#'
#' A functional dataset where the individual functions are given by a Gaussian
#' peak with locations along the \eqn{x}-axis. The variables are as follows: f
#' containing the 29 functions of 101 samples and time which describes the
#' sampling which as been aligned.
#'
#' @format ## `toy_warp`
#' A list which contains the output of the [time_warping()] function applied on
#' the data set `toy_data`.
#'
#' @source ???
"toy_warp"

#' Aligned Simulated two Gaussian Dataset
#'
#' A functional dataset where the individual functions are given by: \eqn{y_i(t)
#' = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3],
#' ~i=1,2,\dots, 21}, where \eqn{z_{i,1}} and \eqn{z_{i,2}} are *i.i.d.* normal
#' with mean one and standard deviation 0.25. Each of these functions is then
#' warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} -
#' 1}) - 3} if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}}
#' (\eqn{gamma_{id}(t) = t}) is the identity warping). The variables are as
#' follows: f containing the 21 functions of 101 samples and time which
#' describes the sampling which has been aligned.
#'
#' @format ## `simu_warp`
#' A list which contains the output of the [time_warping()] function applied on
#' the data set `simu_data`.
#'
#' @source ???
"simu_warp"

#' Aligned Simulated two Gaussian Dataset using Median
#'
#' A functional dataset where the individual functions are given by: \eqn{y_i(t)
#' = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3],
#' ~i=1,2,\dots, 21}, where \eqn{z_{i,1}} and \eqn{z_{i,2}} are *i.i.d.* normal
#' with mean one and standard deviation 0.25. Each of these functions is then
#' warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} -
#' 1}) - 3} if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}}
#' (\eqn{gamma_{id}(t) = t}) is the identity warping). The variables are as
#' follows: f containing the 21 functions of 101 samples and time which
#' describes the sampling which has been aligned.
#'
#' @format ## `simu_warp_median`
#' A list which contains the output of the [time_warping()] function finding the
#' median applied on the data set `simu_data`.
#'
#' @source ???
"simu_warp_median"

#' MPEG7 Curve Dataset
#'
#' Contains the MPEG7 curve data set.
#'
#' @format ## `beta`
#' An array of shape \eqn{2 \times 100 \times 65 \times 20} storing a sample of
#' \eqn{20} curves from \eqn{\mathbb{R}} to \eqn{\mathbb{R}^2} distributed in
#' \eqn{65} different classes, evaluated on a grid of size \eqn{100}.
#'
#' @source ???
"beta"

#' Example Image Data set
#'
#' Contains two simulated images for registration.
#'
#' @format ## `im`
#' A list with two components:
#'
#' - `I1`: A numeric matrix of shape \eqn{64 \times 64} storing the 1st image;
#' - `I2`: A numeric matrix of shape \eqn{64 \times 64} storing the 2nd image.
#'
#' @source ???
"im"
