#ifndef BAYESIAN_H
#define BAYESIAN_H

#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List dpcode(arma::vec q1, arma::vec q1L, arma::vec q2L, int times, int cut);

// [[Rcpp::export]]
Rcpp::List simucode(int iter, int p, arma::vec qt1_5, arma::vec qt2_5, int L,
                    float tau, int times, float kappa, float alpha, float beta,
                    float powera, float dist, float dist_min,
                    arma::vec best_match, arma::vec match, int thin, int cut);

// [[Rcpp::export]]
Rcpp::List itercode(int iter, int n, int m, arma::vec mu_5,
                    arma::mat match_matrix, arma::mat qt_matrix,
                    arma::mat qt_fitted_matrix, int L, float tau, int times,
                    float kappa, float alpha, float beta, float powera,
                    arma::vec best_vec, arma::vec dist_vec,
                    arma::mat best_match_matrix, arma::vec mu_prior,
                    float var_const, arma::vec sumdist, int thin,
                    arma::mat mu_q, arma::mat mu_q_standard, float logmax,
                    int burnin, float AVG);

arma::vec approx(int nd, arma::vec xd, arma::vec yd,int ni, arma::vec xi);
arma::vec findinv(arma::mat warps, int times);
arma::vec R_diff(arma::vec x);

#endif /* BAYESIAN_H */
