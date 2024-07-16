#include <R.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
RcppExport SEXP _fdasrvf_rlbfgs_optim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP calcY(SEXP, SEXP);
RcppExport SEXP check_cross(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP cuL2norm2(SEXP, SEXP);
RcppExport SEXP dpcode(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP DPQ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP DPQ2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP find_grad_2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP find_phistar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP itercode(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP mlogit_warp_grad_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP order_l2norm(SEXP, SEXP);
RcppExport SEXP simucode(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP trapzCpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_fdasrvf_rlbfgs_optim", (DL_FUNC) &_fdasrvf_rlbfgs_optim,  6},
  {"calcY",                 (DL_FUNC) &calcY,                  2},
  {"check_cross",           (DL_FUNC) &check_cross,            4},
  {"cuL2norm2",             (DL_FUNC) &cuL2norm2,              2},
  {"dpcode",                (DL_FUNC) &dpcode,                 5},
  {"DPQ",                   (DL_FUNC) &DPQ,                    8},
  {"DPQ2",                  (DL_FUNC) &DPQ2,                  16},
  {"find_grad_2D",          (DL_FUNC) &find_grad_2D,           6},
  {"find_phistar",          (DL_FUNC) &find_phistar,           7},
  {"itercode",              (DL_FUNC) &itercode,              26},
  {"mlogit_warp_grad_wrap", (DL_FUNC) &mlogit_warp_grad_wrap, 13},
  {"order_l2norm",          (DL_FUNC) &order_l2norm,           2},
  {"simucode",              (DL_FUNC) &simucode,              17},
  {"trapzCpp",              (DL_FUNC) &trapzCpp,               2},
  {NULL, NULL, 0}
};

void attribute_visible R_init_fdasrvf(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
