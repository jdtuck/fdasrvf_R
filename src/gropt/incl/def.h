#ifndef DEF_H
#define DEF_H

#define INITIALBLOCKSIZE 64

//#define MATLAB_MEX_FILE //-------

#ifndef MATLAB_MEX_FILE
	// test examples
//#define TESTEUCFRECHETMEAN
//#define TESTEUCQUADRATIC
//#define TESTPRODUCT
//#define TESTSPHERERAYQUO
//#define TESTSTIEBROCKETT
#define TESTSTIESOFTICA
//#define TESTTESTSPARSEPCA
//#define TESTWEIGHTEDLOWRANK

//#define TESTSIMPLEEXAMPLE
//#define TESTPRODUCTEXAMPLE
#include <stdio.h>
#include <stdlib.h>
// blas and lapack related
#include <cstddef>
#ifdef __cplusplus
    extern "C" {
#endif

extern void dgemm_(
    const char   *transa,
    const char   *transb,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const ptrdiff_t *k,
    const double *alpha,
    const double *a,
    const ptrdiff_t *lda,
    const double *b,
    const ptrdiff_t *ldb,
    const double *beta,
    double *c,
    const ptrdiff_t *ldc
);

extern void dgetrf_(
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    double *a,
    const ptrdiff_t *lda,
    ptrdiff_t *ipiv,
    ptrdiff_t *info
);

extern void dgetrs_(
    const char   *trans,
    const ptrdiff_t *n,
    const ptrdiff_t *nrhs,
    const double *a,
    const ptrdiff_t *lda,
    const ptrdiff_t *ipiv,
    double *b,
    const ptrdiff_t *ldb,
    ptrdiff_t *info
);

extern void dgemv_(
    const char   *trans,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const double *alpha,
    const double *a,
    const ptrdiff_t *lda,
    const double *x,
    const ptrdiff_t *incx,
    const double *beta,
    double *y,
    const ptrdiff_t *incy
);

extern void dcopy_(
    const ptrdiff_t *n,
    const double *dx,
    const ptrdiff_t *incx,
    double *dy,
    const ptrdiff_t *incy
);

extern double ddot_(
    const ptrdiff_t *n,
    const double *dx,
    const ptrdiff_t *incx,
    const double *dy,
    const ptrdiff_t *incy
);

extern void dscal_(
    const ptrdiff_t *n,
    const double *da,
    double *dx,
    const ptrdiff_t *incx
);

extern void daxpy_(
    const ptrdiff_t *n,
    const double *da,
    const double *dx,
    const ptrdiff_t *incx,
    double *dy,
    const ptrdiff_t *incy
);

extern void dger_(
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const double *alpha,
    const double *x,
    const ptrdiff_t *incx,
    const double *y,
    const ptrdiff_t *incy,
    double *a,
    const ptrdiff_t *lda
);

extern void dgeqp3_(
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    double *a,
    const ptrdiff_t *lda,
    ptrdiff_t *jpvt,
    double *tau,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *info
);

extern void dorgqr_(
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const ptrdiff_t *k,
    double *a,
    const ptrdiff_t *lda,
    const double *tau,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *info
);

extern void dormqr_(
    const char   *side,
    const char   *trans,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const ptrdiff_t *k,
    const double *a,
    const ptrdiff_t *lda,
    const double *tau,
    double *c,
    const ptrdiff_t *ldc,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *info
);

extern void dtrsm_(
    const char   *side,
    const char   *uplo,
    const char   *transa,
    const char   *diag,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const double *alpha,
    const double *a,
    const ptrdiff_t *lda,
    double *b,
    const ptrdiff_t *ldb
);

extern void dlarfx_(
    const char   *side,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    const double *v,
    const double *tau,
    double *c,
    const ptrdiff_t *ldc,
    double *work
);

extern void dgesdd_(
    const char   *jobz,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    double *a,
    const ptrdiff_t *lda,
    double *s,
    double *u,
    const ptrdiff_t *ldu,
    double *vt,
    const ptrdiff_t *ldvt,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *iwork,
    ptrdiff_t *info
);

extern void dgesvd_(
    const char   *jobu,
    const char   *jobvt,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    double *a,
    const ptrdiff_t *lda,
    double *s,
    double *u,
    const ptrdiff_t *ldu,
    double *vt,
    const ptrdiff_t *ldvt,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *info
);

extern void dsymv_(
    const char   *uplo,
    const ptrdiff_t *n,
    const double *alpha,
    const double *a,
    const ptrdiff_t *lda,
    const double *x,
    const ptrdiff_t *incx,
    const double *beta,
    double *y,
    const ptrdiff_t *incy
);

extern void dgetri_(
    const ptrdiff_t *n,
    double *a,
    const ptrdiff_t *lda,
    const ptrdiff_t *ipiv,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *info
);

extern void dlapmt_(
    const ptrdiff_t *forwrd,
    const ptrdiff_t *m,
    const ptrdiff_t *n,
    double *x,
    const ptrdiff_t *ldx,
    ptrdiff_t *k
);

extern void dgees_(
    const char   *jobvs,
    const char   *sort,
    ptrdiff_t (*select)(),
    const ptrdiff_t *n,
    double *a,
    const ptrdiff_t *lda,
    ptrdiff_t *sdim,
    double *wr,
    double *wi,
    double *vs,
    const ptrdiff_t *ldvs,
    double *work,
    const ptrdiff_t *lwork,
    ptrdiff_t *bwork,
    ptrdiff_t *info
);

#ifdef __cplusplus
    }   /* extern "C" */
#endif

#endif // end of ifndef MATLAB_MEX_FILE

#ifdef _WIN64
      //define something for Windows (64-bit only)
// Test memory leaking
#ifdef _DEBUG
#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)
//#define CHECKMEMORYDELETED
#else
#define DEBUG_CLIENTBLOCK
#endif

#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

#ifdef _DEBUG
#define new DEBUG_CLIENTBLOCK
#endif

#elif _WIN32
   //define something for Windows (32-bit and 64-bit, this part is common)
#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_IPHONE_SIMULATOR
         // iOS Simulator
    #elif TARGET_OS_IPHONE
        // iOS device
    #elif TARGET_OS_MAC
        // Other kinds of Mac OS
    #else
        // Unsupported platform
    #endif
#elif __linux
    // linux
#elif __unix // all unices not caught above
    // Unix
#elif __posix
    // POSIX
#endif // end of checking platforms
#define integer std::ptrdiff_t

#ifdef MATLAB_MEX_FILE
	#include "mex.h"
    #include "blas.h"
    #include "lapack.h"
    #define integer ptrdiff_t
#define dgemm_ dgemm
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define dgemv_ dgemv
#define dcopy_ dcopy
#define ddot_ ddot
#define dscal_ dscal
#define daxpy_ daxpy
#define dger_ dger
#define dgeqp3_ dgeqp3
#define dorgqr_ dorgqr
#define dormqr_ dormqr
#define dtrsm_ dtrsm
#define dlarfx_ dlarfx
#define dgesdd_ dgesdd
#define dgesvd_ dgesvd
#define dsymv_ dsymv
#define dgetri_ dgetri
#define dgees_ dgees
#endif // end of ifdef MATLAB_MEX_FILE

#include "ForDebug.h"
#include <climits>
#include <limits>
#include "Timer.h"

#ifdef __INTEL_COMPILER
const class {
public:
	template<class T> // convertible to any type
	operator T*(void) const // of null non-member
	{
		return 0;
	} // pointer...
	template<class C, class T> // or any type of null
	operator T C::*(void) const // member pointer...
	{
		return 0;
	}
private:
	void operator&(void) const; // whose address can't be taken
} nullptr = {};
#endif // end of __GNUC__

#include <map>
#include <string>

typedef std::map<std::string, double> PARAMSMAP;

#endif // end of DEF_H
