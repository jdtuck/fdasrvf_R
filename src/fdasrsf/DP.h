#ifndef DP_H
#define DP_H

// pen1: 0 = none, 1 = roughness, 2 = l2gam, 3 = l2psi, 4 = geodesic
// returns 0 on success, -1 if memory could not be allocated
int DP(double *q1, double *q2, int *n1, int *N1, double *lam1, int *pen1,
       int *Disp, double *yy);

#endif /* DP_H */
