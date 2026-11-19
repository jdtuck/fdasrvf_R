/*
 * bayes_warp_geom.cpp  Compiled kernel for the Bayesian HMC alignment
 *
 * Negative log-likelihood and exact gradient (with respect to the basis
 * coefficients c) of
 *
 *   q1(t_k) = (q2 o gam)(t_k) sqrt(gam'(t_k)) + e_k,   e_k ~ N(0, var1)
 *
 * on the uniform grid linspace(0, 1, T). rep = 0: shooting vector
 * (v-transform), rep = 1: clr (h-transform).
 *
 * Batched: c (D,N), q2 (T,N) and q1 (T,1) shared or (T,N) give nll and SSE
 * of length N, grad (D,N), gam and q2w (T,N); each column is an independent
 * warp.
 *
 * This is a port of bayes_warp_geom_c.c from the MATLAB package and
 * implements exactly the same arithmetic.
 */
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

static const double BWG_PI = 3.14159265358979323846;

/* one warp; returns nll, fills grad (if non-NULL), *SSEo, gam, q2w;
 * work must hold 8*T doubles */
static double geom1(const double *c, const double *Phi, arma::uword T,
                    arma::uword D, int rep, const double *q1, const double *q2,
                    double var1, double *grad, double *SSEo, double *gam,
                    double *q2w, double *work)
{
  const double dt = 1.0 / (double)(T - 1);
  double *v = work, *psi = work + T, *apsi = work + 2 * T, *sgn = work + 3 * T;
  double *q2g = work + 4 * T, *q2dg = work + 5 * T, *res = work + 6 * T;
  double *beta = work + 7 * T;

  for (arma::uword k = 0; k < T; k++) v[k] = 0.0;
  for (arma::uword j = 0; j < D; j++) {
    const double cj = c[j];
    const double *col = Phi + j * T;
    for (arma::uword k = 0; k < T; k++) v[k] += col[k] * cj;
  }

  double r = 0.0, s = 1.0, ds = 0.0, sr = 0.0;
  if (rep == 0) {
    for (arma::uword k = 0; k < T; k++) {
      double wk = (k == 0 || k == T - 1) ? dt / 2 : dt;
      r += wk * v[k] * v[k];
    }
    r = std::sqrt(r);
    if (r >= BWG_PI) {
      for (arma::uword k = 0; k < T; k++) {
        gam[k] = k * dt;
        q2w[k] = q2[k];
      }
      if (grad) for (arma::uword j = 0; j < D; j++) grad[j] = 0.0;
      *SSEo = R_PosInf;
      return R_PosInf;
    }
    double cr = 1.0;
    if (r >= 1e-10) {
      s = std::sin(r) / r;
      ds = (r * std::cos(r) - std::sin(r)) / (r * r);
      sr = std::sin(r);
      cr = std::cos(r);
    }
    for (arma::uword k = 0; k < T; k++) {
      psi[k] = cr + s * v[k];
      sgn[k] = (psi[k] < 0) ? -1.0 : 1.0;
      apsi[k] = std::fabs(psi[k]);
    }
  } else {
    double vmax = v[0];
    for (arma::uword k = 1; k < T; k++) if (v[k] > vmax) vmax = v[k];
    for (arma::uword k = 0; k < T; k++) {
      psi[k] = std::exp((v[k] - vmax) / 2);
      sgn[k] = 1.0;
      apsi[k] = psi[k];
    }
  }

  gam[0] = 0.0;
  for (arma::uword k = 1; k < T; k++)
    gam[k] = gam[k - 1] + dt / 2 * (psi[k - 1] * psi[k - 1] + psi[k] * psi[k]);
  const double GT = gam[T - 1];
  for (arma::uword k = 0; k < T; k++) gam[k] /= GT;
  gam[T - 1] = 1.0;

  const double a = 1.0 / std::sqrt(GT);
  double SSE = 0.0, rq = 0.0;
  for (arma::uword k = 0; k < T; k++) {
    double pos = gam[k] * (T - 1);
    if (pos < 0) pos = 0;
    if (pos > (double)(T - 1)) pos = (double)(T - 1);
    arma::uword idx = (arma::uword)std::floor(pos);
    if (idx > T - 2) idx = T - 2;
    double fr = pos - (double)idx;
    q2g[k] = q2[idx] * (1 - fr) + q2[idx + 1] * fr;
    q2dg[k] = (q2[idx + 1] - q2[idx]) * (T - 1);
    q2w[k] = q2g[k] * apsi[k] * a;
    res[k] = q1[k] - q2w[k];
    SSE += res[k] * res[k];
    rq += res[k] * q2w[k];
  }
  *SSEo = SSE;

  if (grad) {
    const double kappa = rq / (2 * GT);
    double Agam = 0.0;
    for (arma::uword k = 0; k < T; k++) {
      double A = res[k] * a * q2dg[k] * apsi[k];
      Agam += A * gam[k];
      q2dg[k] = A;
    }
    double Ssuf = 0.0;
    for (arma::uword kk = T; kk-- > 0;) {
      double Sk1 = Ssuf;
      Ssuf += q2dg[kk];
      double coefP = 0.0;
      if (kk < T - 1) coefP += dt / 2 * Sk1;
      if (kk > 0) coefP += dt / 2 * Ssuf;
      double wk = (kk == 0 || kk == T - 1) ? dt / 2 : dt;
      double alpha = coefP / GT - (Agam / GT + kappa) * wk;
      double Bv = res[kk] * a * q2g[kk] * sgn[kk];
      beta[kk] = 2 * alpha * psi[kk] + Bv;
    }
    if (rep == 0) {
      double sb = 0.0, sbv = 0.0;
      for (arma::uword k = 0; k < T; k++) {
        sb += beta[k];
        sbv += beta[k] * v[k];
      }
      double fac = (r >= 1e-10) ? (-sr * sb + ds * sbv) / r : 0.0;
      for (arma::uword j = 0; j < D; j++) {
        const double *col = Phi + j * T;
        double pb = 0.0, pwv = 0.0;
        for (arma::uword k = 0; k < T; k++) {
          double wk = (k == 0 || k == T - 1) ? dt / 2 : dt;
          pb += col[k] * beta[k];
          pwv += col[k] * wk * v[k];
        }
        grad[j] = -(s * pb + fac * pwv) / var1;
      }
    } else {
      for (arma::uword j = 0; j < D; j++) {
        const double *col = Phi + j * T;
        double pb = 0.0;
        for (arma::uword k = 0; k < T; k++) pb += col[k] * beta[k] * psi[k];
        grad[j] = -0.5 * pb / var1;
      }
    }
  }
  return SSE / (2 * var1) + 0.5 * (double)T * std::log(var1);
}

//' Negative log-likelihood and gradient of a Bayesian warp (internal)
//'
//' @param c matrix of basis coefficients, one warping function per column
//' @param Phi basis matrix on the uniform grid
//' @param rep 0 for the shooting-vector transform, 1 for the clr transform
//' @param q1 target SRVF, either one column or one per column of `c`
//' @param q2 SRVFs to be warped, one per column of `c`
//' @param var1 noise variance
//' @param do_grad whether to compute the gradient
//'
//' @return A list with the negative log-likelihoods, gradients, SSEs,
//'   warping functions and warped SRVFs.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List bayes_warp_geom_c(const arma::mat& c, const arma::mat& Phi, int rep,
                             const arma::mat& q1, const arma::mat& q2,
                             double var1, bool do_grad = true)
{
  const arma::uword T = Phi.n_rows;
  const arma::uword D = Phi.n_cols;
  if (T < 2) Rcpp::stop("need at least two grid points");
  if (c.n_rows != D) Rcpp::stop("nrow(c) must equal ncol(Phi)");
  const arma::uword N = c.n_cols;
  if (q2.n_rows != T || q2.n_cols != N)
    Rcpp::stop("q2 must be (nrow(Phi), ncol(c))");
  bool q1_shared;
  if (q1.n_rows == T && q1.n_cols == 1)
    q1_shared = true;
  else if (q1.n_rows == T && q1.n_cols == N)
    q1_shared = false;
  else
    Rcpp::stop("q1 must be (nrow(Phi), 1) or the size of q2");

  arma::vec nll(N), sse(N);
  arma::mat grad(D, N, arma::fill::zeros);
  arma::mat gam(T, N), q2w(T, N);
  std::vector<double> work(8 * T);

  for (arma::uword i = 0; i < N; i++) {
    nll[i] = geom1(c.colptr(i), Phi.memptr(), T, D, rep,
                   q1_shared ? q1.memptr() : q1.colptr(i), q2.colptr(i), var1,
                   do_grad ? grad.colptr(i) : NULL, &sse[i], gam.colptr(i),
                   q2w.colptr(i), work.data());
  }

  return Rcpp::List::create(Rcpp::Named("nll") = nll,
                            Rcpp::Named("grad") = grad,
                            Rcpp::Named("SSE") = sse,
                            Rcpp::Named("gam") = gam,
                            Rcpp::Named("q2w") = q2w);
}
