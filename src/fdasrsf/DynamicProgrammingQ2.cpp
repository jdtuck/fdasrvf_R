#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>
#include "dp_grid.h"
#include "DynamicProgrammingQ2.h"
using namespace Rcpp;


RcppExport SEXP DPQ2(SEXP Q1, SEXP T1, SEXP Q2, SEXP T2, SEXP m1, SEXP n1, SEXP n2, SEXP tv1, SEXP tv2, SEXP n1v, SEXP n2v, SEXP G, SEXP T, SEXP size, SEXP lam1){

  NumericVector Q1i(Q1);
  NumericVector Q2i(Q2);
  NumericVector T1i(T1);
  NumericVector T2i(T2);
  NumericVector tv1i(tv1);
  NumericVector tv2i(tv2);
  NumericVector GG(G);
  NumericVector TT(T);

  double * _Q1i = &Q1i[0];
  double * _Q2i = &Q2i[0];
  double * _T1i = &T1i[0];
  double * _T2i = &T2i[0];
  double * _tv1i = &tv1i[0];
  double * _tv2i = &tv2i[0];
  double * _GG = &GG[0];
  double * _TT = &TT[0];

  int _m1 = as<int>(m1);
  int _n1 = as<int>(n1);
  int _n2 = as<int>(n2);
  int _n1v = as<int>(n1v);
  int _n2v = as<int>(n2v);
  int _size = as<int>(size);
  double _lam1 = as<double>(lam1);

  DynamicProgrammingQ2(_Q1i, _T1i, _Q2i, _T2i, &_m1, &_n1, &_n2, _tv1i, _tv2i, &_n1v, &_n2v, _GG, _TT, &_size, &_lam1);

  List ret; ret["G"] = wrap(GG); ret["T"] = wrap(TT); ret["size"] = wrap(_size);
  return(ret);
}

void DynamicProgrammingQ2(double *Q1, double *T1, double *Q2, double *T2, const int *m1, const int *n1, const int *n2, double *tv1, double *tv2, const int *n1v, const int *n2v, double *G, double *T, int *size, const double *lam1){
  int *idxv1 = 0;
  int *idxv2 = 0;
  double *E = 0; /* E[ntv1*j+i] = cost of best path to (tv1[i],tv2[j]) */
  int *P = 0; /* P[ntv1*j+i] = predecessor of (tv1[i],tv2[j]) along best path */

  idxv1=(int*)malloc((*n1v)*sizeof(int));
  idxv2=(int*)malloc((*n2v)*sizeof(int));
  E=(double*)malloc((*n1v)*(*n2v)*sizeof(double));
  P=(int*)calloc((*n1v)*(*n2v),sizeof(int));

  /* dp_costs() needs indexes for gridpoints precomputed */
  dp_all_indexes( T1, *n1, tv1, *n1v, idxv1 );
  dp_all_indexes( T2, *n2, tv2, *n2v, idxv2 );

  /* Compute cost of best path from (0,0) to every other grid point */
  dp_costs( Q1, T1, *n1, Q2, T2, *n2,
    *m1, tv1, idxv1, *n1v, tv2, idxv2, *n2v, E, P, *lam1 );

  /* Reconstruct best path from (0,0) to (1,1) */
  *size = dp_build_gamma( P, tv1, *n1v, tv2, *n2v, G, T );

  // free allocated memory
  free(idxv1); free(idxv2); free(E); free(P);
}
