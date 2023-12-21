#include "fdasrsf/rbfgs.h"
#include <RcppArmadillo.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

RcppExport SEXP _fdasrvf_rlbfgs_optim(SEXP q1SEXP, SEXP q2SEXP, SEXP timeSEXP, SEXP maxiterSEXP, SEXP lamSEXP, SEXP penaltySEXP) {
  vec q1 = as<vec>(q1SEXP);
  vec q2 = as<vec>(q2SEXP);
  vec time = as<vec>(timeSEXP);
  int maxiter = as<int>(maxiterSEXP);
  int len = q1.size();
  double lam = as<double>(lamSEXP);
  int penalty = as<int>(penaltySEXP);
  vec gam(len);
  gam = rlbfgs_optim(q1, q2, time, maxiter, lam, penalty);
  return wrap(gam);
}
