#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

int main() {
    uword T = 101;
    vec time = arma::linspace(0, 2*M_PI, T);
    vec q1 = sin(time);
    vec q2 = cos(time);
    vec time1 = arma::linspace(0, 1, T);

    rlbfgs myObj(q1, q2, time1);
    myObj.solve();

    myObj.gammaOpt.print();

    return 0;
}
