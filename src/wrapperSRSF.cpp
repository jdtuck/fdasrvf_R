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
                Rcpp::NumericVector G,
                Rcpp::NumericVector T,
                int size,
                double lam1,
                int nbhd_dim)
{
  DynamicProgrammingQ2(Q1.begin(), T1.begin(), Q2.begin(), T2.begin(), &m1, &n1,
                       &n2, tv1.begin(), tv2.begin(), &n1v, &n2v, G.begin(),
                       T.begin(), &size, &lam1, &nbhd_dim);

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
  Rcpp::NumericVector out(N1);
  DP(Q1.begin(), Q2.begin(), &n1, &N1, &lam1, &pen1, &Disp, out.begin());
  return(out);
}

// [[Rcpp::export]]
arma::vec rlbfgs(arma::vec q1, arma::vec q2, arma::vec time,
                                int maxiter, double lam, int penalty) {
  arma::vec gam = rlbfgs_optim(q1, q2, time, maxiter, lam, penalty);
  return(gam);
}
