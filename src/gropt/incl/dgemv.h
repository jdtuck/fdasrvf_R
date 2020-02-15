#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 
#define dgemv1_  dgemv

int dgemv(char *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy, int trans_len);

#ifdef __cplusplus
}
#endif