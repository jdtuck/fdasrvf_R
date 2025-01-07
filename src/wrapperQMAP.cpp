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
Rcpp::NumericVector find_phistar(Rcpp::NumericVector w,
                                 Rcpp::NumericVector q,
                                 Rcpp::NumericVector b,
                                 int n, int t, int d, int K) {
  findphistar(w.begin(), q.begin(), b.begin(), n, t, d, K);
  return(w);
}
