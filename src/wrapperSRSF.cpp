#include "fdasrsf/mlogit_warp_grad.h"
#include "fdasrsf/DynamicProgrammingQ2.h"
#include "fdasrsf/DP.h"
#include "fdasrsf/rbfgs.h"
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::NumericVector mlogit_warp_grad_wrap(int m1, int m2,
                                          Rcpp::NumericVector alpha,
                                          Rcpp::NumericVector beta,
                                          Rcpp::NumericVector ti,
                                          Rcpp::NumericVector gami,
                                          Rcpp::NumericVector q,
                                          Rcpp::IntegerVector y, int max_itri,
                                          double toli, double deltai,
                                          int displayi) {
  Rcpp::NumericVector gamout(m1);
  mlogit_warp_grad(&m1, &m2, alpha.begin(), beta.begin(), ti.begin(),
                   gami.begin(), q.begin(), y.begin(), &max_itri, &toli,
                   &deltai, &displayi, gamout.begin());
  return gamout;
}

// [[Rcpp::export]]
Rcpp::List DPQ2(Rcpp::NumericVector Q1,
                Rcpp::NumericVector T1,
                Rcpp::NumericVector Q2,
                Rcpp::NumericVector T2,
                int m1,
                int n1,
                int n2,
                Rcpp::NumericVector tv1,
                Rcpp::NumericVector tv2,
                int n1v,
                int n2v,
                double lam1,
                int nbhd_dim,
                int pen = 1)
{
  // 0 = no penalty, 1 = roughness, 2 = l2gam, 3 = l2psi, 4 = geodesic
  if (pen < 0 || pen > 4)
    Rcpp::stop("pen must be one of 0, 1, 2, 3 or 4.");
  if (nbhd_dim < 1 || nbhd_dim > 65535)
    Rcpp::stop("nbhd_dim must be between 1 and 65535.");
  if (m1 < 1 || n1 < 2 || n2 < 2 || n1v < 2 || n2v < 2)
    Rcpp::stop("m1 must be positive and n1, n2, n1v and n2v at least two.");
  if (T1.size() < n1 || T2.size() < n2 || tv1.size() < n1v || tv2.size() < n2v)
    Rcpp::stop("T1, T2, tv1 and tv2 are shorter than n1, n2, n1v and n2v.");
  // dp_edge_weight() interpolates Q linearly, reading columns 0 .. n-1
  if (Q1.size() < (R_xlen_t)m1 * n1 || Q2.size() < (R_xlen_t)m1 * n2)
    Rcpp::stop("Q1 and Q2 must have at least m1*n1 and m1*n2 elements.");

  // dp_build_gamma() needs room for max(n1v, n2v) points
  int Gsize = n1v > n2v ? n1v : n2v;
  Rcpp::NumericVector G(Gsize);
  Rcpp::NumericVector T(Gsize);
  int size = 0;
  if (DynamicProgrammingQ2(Q1.begin(), T1.begin(), Q2.begin(), T2.begin(), &m1,
                           &n1, &n2, tv1.begin(), tv2.begin(), &n1v, &n2v,
                           G.begin(), T.begin(), &size, &lam1, &nbhd_dim,
                           &pen) != 0)
    Rcpp::stop("DynamicProgrammingQ2: out of memory.");

  Rcpp::List ret;
  ret["G"] = G;
  ret["T"] = T;
  ret["size"] = size;
  return(ret);
}

// [[Rcpp::export]]
Rcpp::NumericVector DPQ(Rcpp::NumericVector Q1,
                        Rcpp::NumericVector Q2,
                        int n1,
                        int N1,
                        double lam1,
                        int pen1,
                        int Disp,
                        int scl = 3) {
  // 0 = no penalty, 1 = roughness, 2 = l2gam, 3 = l2psi, 4 = geodesic
  if (pen1 < 0 || pen1 > 4)
    Rcpp::stop("pen1 must be one of 0, 1, 2, 3 or 4.");
  // q1 and q2 are resampled scl times finer to evaluate the cost
  if (scl < 1 || scl > 64)
    Rcpp::stop("scl must be an integer between 1 and 64.");
  if (n1 < 1 || N1 < 2)
    Rcpp::stop("n1 must be positive and N1 at least two.");
  if (Q1.size() < (R_xlen_t)n1 * N1 || Q2.size() < (R_xlen_t)n1 * N1)
    Rcpp::stop("Q1 and Q2 must each have at least n1*N1 elements.");

  Rcpp::NumericVector out(N1);
  if (DP(Q1.begin(), Q2.begin(), &n1, &N1, &lam1, &pen1, &scl, &Disp, out.begin()) != 0)
    Rcpp::stop("DP: out of memory.");
  return(out);
}

// gam0, when not NULL, is the warping to start from instead of the identity.
// info = [iterations, stop reason, cost, gradient norm]; the stop reason is
// 1 for no acceptable step, 2 for the gradient norm, 3 for maxiter and 4
// for the relative change of the cost
// [[Rcpp::export]]
Rcpp::List rlbfgs(arma::vec q1, arma::vec q2, arma::vec time,
                  int maxiter, double lam, int penalty,
                  Rcpp::Nullable<Rcpp::NumericVector> gam0 = R_NilValue) {
  if (q2.n_elem != q1.n_elem || time.n_elem != q1.n_elem)
    Rcpp::stop("q1, q2 and time must have the same number of elements.");
  arma::vec g0;
  if (gam0.isNotNull()) {
    g0 = Rcpp::as<arma::vec>(gam0.get());
    if (!g0.is_empty() && g0.n_elem != time.n_elem)
      Rcpp::stop("gam0 must have as many elements as time.");
  }
  arma::vec info;
  arma::vec gam = rlbfgs_optim(q1, q2, time, maxiter, lam, penalty, g0, info);
  return Rcpp::List::create(Rcpp::Named("gam") = gam,
                            Rcpp::Named("info") = info);
}
