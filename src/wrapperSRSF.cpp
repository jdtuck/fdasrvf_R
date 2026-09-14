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
  // dp_edge_weight() reads Q columns 0 .. n-2
  if (Q1.size() < (R_xlen_t)m1 * (n1 - 1) || Q2.size() < (R_xlen_t)m1 * (n2 - 1))
    Rcpp::stop("Q1 and Q2 must have at least m1*(n1-1) and m1*(n2-1) elements.");

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
                        int Disp) {
  // 0 = no penalty, 1 = roughness, 2 = l2gam, 3 = l2psi, 4 = geodesic
  if (pen1 < 0 || pen1 > 4)
    Rcpp::stop("pen1 must be one of 0, 1, 2, 3 or 4.");
  if (n1 < 1 || N1 < 2)
    Rcpp::stop("n1 must be positive and N1 at least two.");
  if (Q1.size() < (R_xlen_t)n1 * N1 || Q2.size() < (R_xlen_t)n1 * N1)
    Rcpp::stop("Q1 and Q2 must each have at least n1*N1 elements.");

  Rcpp::NumericVector out(N1);
  if (DP(Q1.begin(), Q2.begin(), &n1, &N1, &lam1, &pen1, &Disp, out.begin()) != 0)
    Rcpp::stop("DP: out of memory.");
  return(out);
}

// [[Rcpp::export]]
arma::vec rlbfgs(arma::vec q1, arma::vec q2, arma::vec time,
                                int maxiter, double lam, int penalty) {
  arma::vec gam = rlbfgs_optim(q1, q2, time, maxiter, lam, penalty);
  return(gam);
}
