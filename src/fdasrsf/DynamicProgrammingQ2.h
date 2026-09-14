#ifndef DYNAMICPROGRAMMINGQ2_H
#define DYNAMICPROGRAMMINGQ2_H

// G and T must have room for max(n1v, n2v) points.
// pen1: 0 = none, 1 = roughness, 2 = l2gam, 3 = l2psi, 4 = geodesic
// returns 0 on success, -1 if memory could not be allocated
int DynamicProgrammingQ2(double *Q1, double *T1, double *Q2, double *T2,
                         const int *m1, const int *n1, const int *n2,
                         double *tv1, double *tv2, const int *n1v,
                         const int *n2v, double *G, double *T, int *size,
                         const double *lam1, int *nbhd_dim1, const int *pen1);

#endif /* DYNAMICPROGRAMMINGQ2_H */
