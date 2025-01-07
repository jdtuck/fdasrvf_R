#include "pair_align_functions_expomap.h"

arma::vec calcY(double area, arma::vec gy) {
  int len = gy.size();
  double sarea = std::sin(area);
  double carea = std::cos(area);
  arma::vec output(len);
  if (area != 0) {
    for (unsigned int i = 0; i < len; i++) {
      output[i] = carea + sarea / area * gy[i];
    }
    return output;
  } else {
    for (unsigned int i = 0; i < len; i++) {
      output[i] = 1.0;
    }
    return output;
  }
}

Rcpp::NumericVector cuL2norm2(arma::vec x, arma::vec y) {
  int n = x.size();
  // get ordering of x
  // borrowed from Hadley Wickhams github.com/hadley/adv-r/blob/master/extras/cpp/order.cpp
  std::vector<std::pair<double, int> > vals;
  vals.reserve(n);
  for (unsigned int i = 0; i < n; i++) {
    vals.push_back(std::make_pair(x[i], i));
  }
  std::sort(vals.begin(), vals.end());
  int sInd;
  Rcpp::NumericVector xSort(n);
  Rcpp::NumericVector ySortSq (n);
  for (unsigned int i = 0; i < n; i++) {
    sInd = vals[i].second;
    xSort[i] = x[sInd];
    ySortSq[i] = std::pow(y[sInd], 2.0);
  }
  // Get trapezoid areas
  Rcpp::NumericVector prod_xy2(n - 1);
  for (unsigned int i = 0; i < (n - 1); i++) {
    prod_xy2[i] = (xSort[i+1] - xSort[i]) * (ySortSq[i+1] + ySortSq[i]);
  }
  // cumulative sum
  Rcpp::NumericVector cusum(n - 1);
  for (unsigned int i = 0; i < (n - 1); i++) {
    cusum[i] = 0;
    for (unsigned int j = 0; j <= i; j++) {
      cusum[i] += prod_xy2[j];
    }
  }
  return cusum / 2.0;
}

double trapzCpp(arma::vec x, arma::vec y) {
  int n = x.size();
  double area2 = 0.0;
  for (int i = 0; i < (n-1); i++) {
    area2 += (x[i+1] - x[i]) * (y[i+1] + y[i]);
  }
  return area2 / 2.0;
}

double order_l2norm(arma::vec x, arma::vec y) {
  int n = x.size();
  // get ordering of x
  // borrowed from Hadley Wickhams adv-r/extras/cpp/order.cpp
  std::vector<std::pair<double, int> > vals;
  vals.reserve(n);
  for (unsigned int i = 0; i < n; i++) {
    vals.push_back(std::make_pair(x[i], i));
  }
  std::sort(vals.begin(), vals.end());
  int sInd;
  Rcpp::NumericVector xSort(n);
  Rcpp::NumericVector ySortSq (n);
  for (unsigned int i = 0; i < n; i++) {
    sInd = vals[i].second;
    xSort[i] = x[sInd];
    ySortSq[i] = std::pow(y[sInd], 2.0);
  }
  // loop through sorted inds, square y, trapz
  double area2 = 0.0;
  for (unsigned int i = 0; i < (n - 1); i++) {
    area2 += (xSort[i+1] - xSort[i]) * (ySortSq[i+1] + ySortSq[i]);
  }
  return std::sqrt(area2 / 2.0);
}
