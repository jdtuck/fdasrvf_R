\name{horizFPCA}
\alias{horizFPCA}
\title{Horizontal Functional Principal Component Analysis}
\usage{
  horizFPCA(gam, no, showplot = TRUE)
}
\arguments{
  \item{gam}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} of
  warping functions}

  \item{no}{number of prinicpal components to extract}

  \item{showplot}{show plots of prinipal directions
  (default = T)}
}
\value{
  Returns a list containing \item{gam_pca}{warping
  functions principal directions} \item{psi_pca}{srvf
  principal directions} \item{latent}{latent values}
  \item{U}{eigenvectors}
}
\description{
  This function calculates vertical functional principal
  component analysis on aligned data
}
\examples{
data("simu_data")
gam = time_warping(simu_data$f,simu_data$time)$gam
out = SqrtMean(gam)
hfpca = horizFPCA(out$vec,out$mu,no = 3)
}
\references{
  Tucker, J. D., Wu, W., Srivastava, A., May 2012.
  Generative Models for Function Data using Phase and
  Amplitude Separation, submitted to Computational
  Statistics and Data Analysis.
}
\keyword{alignment}
\keyword{srvf}

