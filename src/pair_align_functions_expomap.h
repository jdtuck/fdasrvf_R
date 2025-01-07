#include <RcppArmadillo.h>

// calculate exponential map in f.exp1()
// [[Rcpp::export]]
arma::vec calcY(double area, arma::vec gy);

// cumulative sum of squares for f.SSEg.pw()
// [[Rcpp::export]]
Rcpp::NumericVector cuL2norm2(arma::vec x, arma::vec y);

// simple trapezoidal numerical integration
// [[Rcpp::export]]
double trapzCpp(arma::vec x, arma::vec y);

// order vectors and calculate l2 norm, for f.L2norm()
// [[Rcpp::export]]
double order_l2norm(arma::vec x, arma::vec y);
