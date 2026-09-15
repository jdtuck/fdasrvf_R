#include "fdaqmap/incl/UnitSquareImage.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List find_grad_2D(Rcpp::NumericVector dfdu,
                        Rcpp::NumericVector dfdv,
                        Rcpp::NumericVector f,
                        int n, int t, int d) {
  findgrad2D(dfdu.begin(), dfdv.begin(), f.begin(), n, t, d);
  Rcpp::List ret;
  ret["dfdu"] = dfdu;
  ret["dfdv"] = dfdv;
  return(ret);
}

// [[Rcpp::export]]
int check_cross(Rcpp::NumericVector f, int n, int t, int D) {
  return check_crossing(f.begin(), n, t, D);
}

// [[Rcpp::export]]
Rcpp::NumericVector interp_surf(Rcpp::NumericVector F,
                                Rcpp::NumericVector u,
                                Rcpp::NumericVector v,
                                int m, int n, int d) {
  int P = u.size();
  if (v.size() != P)
    Rcpp::stop("u and v must have the same length");
  if (m < 2 || n < 2 || d < 1 || F.size() != m*n*d)
    Rcpp::stop("F must be an m x n x d array with m, n >= 2");
  Rcpp::NumericVector Fnew(P*d);
  Interp_Surf(Fnew.begin(), F.begin(), u.begin(), v.begin(), m, n, d, P);
  return(Fnew);
}

// [[Rcpp::export]]
Rcpp::NumericVector find_phistar(Rcpp::NumericVector w,
                                 Rcpp::NumericVector q,
                                 Rcpp::NumericVector b,
                                 int n, int t, int d, int K) {
  findphistar(w.begin(), q.begin(), b.begin(), n, t, d, K);
  return(w);
}
