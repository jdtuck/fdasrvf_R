#include <RcppArmadillo.h>
#include <Rcpp.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// calculate exponential map in f.exp1()
RcppExport SEXP calcY(SEXP R_area, SEXP R_gy) {
  vec gy = as<vec>(R_gy);
  int len = gy.size();
  double area = as<double>(R_area);
  double sarea = sin(area);
  double carea = cos(area);
  vec output(len);
  if (area != 0) {
    for (int i = 0; i < len; i++) {
      output[i] = carea + sarea / area * gy[i];
    }
    return wrap(output);
  } else {
    for (int i = 0; i < len; i++) {
      output[i] = 1.0;
    }
    return wrap(output);
  }
}


// cumulative sum of squares for f.SSEg.pw()
//RcppExport NumericVector cuL2norm2(NumericVector x, NumericVector y) {
RcppExport SEXP cuL2norm2(SEXP R_x, SEXP R_y) {
  vec x = as<vec>(R_x);
  vec y = as<vec>(R_y);
  int n = x.size();

  // get ordering of x
  // borrowed from Hadley Wickhams github.com/hadley/adv-r/blob/master/extras/cpp/order.cpp
  
  vector<pair<double, int> > vals;
  vals.reserve(n);
  for(int i = 0; i < n; i++) {
    vals.push_back(make_pair(x[i], i));
  }
  /*
  std::vector<std::pair <double, int> > vals(n);
  for(int i = 0; i < n; i++) {
    vals[i] = std::make_pair<double, int>(x[i], i);
  }
  */

  std::sort(vals.begin(), vals.end());
  int sInd;
  NumericVector xSort(n);
  NumericVector ySortSq (n);
  for(int i = 0; i < n; i++) {
    sInd = vals[i].second;
    xSort[i] = x[sInd];
    ySortSq[i] = pow(y[sInd], 2.0);
  }

  // Get trapezoid areas
  NumericVector prod_xy2(n-1);
  for (int i = 0; i < (n-1); i++) {
    prod_xy2[i] = (xSort[i+1] - xSort[i]) * (ySortSq[i+1] + ySortSq[i]);
  }

  // cumulative sum
  NumericVector cusum(n-1);
  for (int i = 0; i < (n-1); i++) {
    cusum[i] = 0;
    for (int j = 0; j <= i; j++) {
      cusum[i] += prod_xy2[j];
    }
  }
  return wrap(cusum / 2.0);
}

// simple trapezoidal numerical integration
//RcppExport double trapzCpp(NumericVector x, NumericVector y) {
RcppExport SEXP trapzCpp(SEXP R_x, SEXP R_y) {
  vec x = as<vec>(R_x);
  vec y = as<vec>(R_y);
  int n = x.size();
  double area2 = 0.0;
  for (int i = 0; i < (n-1); i++) {
    area2 += (x[i+1] - x[i]) * (y[i+1] + y[i]);
  }
  return wrap(area2/2.0);
}

// order vectors and calculate l2 norm, for f.L2norm()
//RcppExport double order_l2norm(NumericVector x, NumericVector y) {
RcppExport SEXP order_l2norm(SEXP R_x, SEXP R_y) {
  vec x = as<vec>(R_x);
  vec y = as<vec>(R_y);
  int n = x.size();
  // get ordering of x
  // borrowed from Hadley Wickhams adv-r/extras/cpp/order.cpp

  vector<pair<double, int> > vals;
  vals.reserve(n);
  for(int i = 0; i < n; i++) {
    vals.push_back(make_pair(x[i], i));
  }

  /*
  std::vector<std::pair <double, int> > vals(n);
  for(int i = 0; i < n; i++) {
    vals[i] = std::make_pair<double, int>(x[i], i);
  }
  */

  std::sort(vals.begin(), vals.end());
  int sInd;
  NumericVector xSort(n);
  NumericVector ySortSq (n);
  for(int i = 0; i < n; i++) {
    sInd = vals[i].second;
    xSort[i] = x[sInd];
    ySortSq[i] = pow(y[sInd], 2.0);
  }

  // loop through sorted inds, square y, trapz
  double area2 = 0.0;
  for (int i = 0; i<(n-1); i++) {
    area2 += (xSort[i+1] - xSort[i]) * (ySortSq[i+1] + ySortSq[i]);
  }
  return wrap(sqrt(area2 / 2.0));
}
