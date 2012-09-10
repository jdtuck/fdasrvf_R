#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "dp_grid.h"

void DynamicProgrammingQ2(double *Q1, double *T1, double *Q2, double *T2, int *m1, int *n1, int *n2, double *tv1, double *tv2, int *n1v, int *n2v, double *G, double *T, int *size){
  int nsamps1;
  int nsamps2;
  int *idxv1 = 0;
  int *idxv2 = 0;
  int ntv1;
  int ntv2;
  int dim = 0;
  double *E = 0; /* E[ntv1*j+i] = cost of best path to (tv1[i],tv2[j]) */
  int *P = 0; /* P[ntv1*j+i] = predecessor of (tv1[i],tv2[j]) along best path */
  double m, rootm;
  int sr, sc; /* source row and column index */
  int tr, tc; /* target row and column index */
  int Galloc_size;
  double pres = 0;

  dim = *m1;
  nsamps1 = *n1; /* = columns(T1) = columns(Q1)+1 */
  nsamps2 = *n2; /* = columns(T2) = columns(Q2)+1 */
  ntv1 = *n1v;
  ntv2 = *n2v;
  Galloc_size = ntv1>ntv2 ? ntv1 : ntv2;

  idxv1=(int*)malloc(ntv1*sizeof(int));
  idxv2=(int*)malloc(ntv2*sizeof(int));
  E=(double*)malloc(ntv1*ntv2*sizeof(double));
  P=(int*)calloc(ntv1*ntv2,sizeof(int));

  /* dp_costs() needs indexes for gridpoints precomputed */
  dp_all_indexes( T1, nsamps1, tv1, ntv1, idxv1 );
  dp_all_indexes( T2, nsamps2, tv2, ntv2, idxv2 );

  // /* Compute cost of best path from (0,0) to every other grid point */
  pres = dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2, 
    dim, tv1, idxv1, ntv1, tv2, idxv2, ntv2, E, P );

  // /* Reconstruct best path from (0,0) to (1,1) */
  *size = dp_build_gamma( P, tv1, ntv1, tv2, ntv2, G, T );

}
