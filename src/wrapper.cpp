
// #include "gropt/incl/ElasticCurvesReparam.h"
#include "fdasrsf/DynamicProgrammingQ2.h"
#include "fdasrsf/mlogit_warp_grad.h"
#include "fdasrsf/DP.h"
#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP mlogit_warp_grad_wrap(SEXP m1, SEXP m2, SEXP alpha, SEXP beta, SEXP ti, SEXP gami, SEXP q, SEXP y, SEXP max_itri, SEXP toli, SEXP deltai, SEXP displayi, SEXP gamout){
  NumericVector alphai(alpha);
  NumericVector betai(beta);
  NumericVector tii(ti);
  NumericVector gamii(gami);
  NumericVector qi(q);
  IntegerVector yi(y);
  NumericVector gamouti(gamout);

  double * _alphai = &alphai[0];
  double * _betai = &betai[0];
  double * _tii = &tii[0];
  double * _gamii = &gamii[0];
  double * _qi = &qi[0];
  int * _yi = &yi[0];
  double * _gamouti = &gamouti[0];

  int _m1 = as<int>(m1);
  int _m2 = as<int>(m2);
  int _max_itri = as<int>(max_itri);
  double _toli = as<double>(toli);
  double _deltai = as<double>(deltai);
  int _displayi = as<int>(displayi);

  mlogit_warp_grad(&_m1, &_m2, _alphai, _betai, _tii, _gamii, _qi, _yi, &_max_itri, &_toli, &_deltai, &_displayi, _gamouti);

  return(gamouti);
}

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

RcppExport SEXP DPQ(SEXP Q1, SEXP Q2, SEXP n1, SEXP N1, SEXP lam1, SEXP Disp, SEXP yy){

  NumericVector Q1i(Q1);
  NumericVector Q2i(Q2);
  NumericVector yyi(yy);

  double * _Q1i = &Q1i[0];
  double * _Q2i = &Q2i[0];
  double * _yyi = &yyi[0];

  int _n1 = as<int>(n1);
  int _N1 = as<int>(N1);
  int _Disp = as<int>(Disp);
  double _lam1 = as<double>(lam1);

  DP(_Q1i, _Q2i, &_n1, &_N1, &_lam1, &_Disp, _yyi);

  return(yyi);
}

// RcppExport SEXP opt_reparam(SEXP C1, SEXP C2, SEXP n, SEXP d, SEXP w,
//                             SEXP onlyDP, SEXP rotated, SEXP isclosed,
//                             SEXP skipm, SEXP autoselectC, SEXP opt, SEXP swap,
//                             SEXP fopts, SEXP comtime){
// 
//     NumericVector C1i(C1);
//     NumericVector C2i(C2);
//     NumericVector opti(opt);
//     NumericVector foptsi(fopts);
//     NumericVector comtimei(comtime);
// 
//     double *_C1i = &C1i[0];
//     double *_C2i = &C2i[0];
//     double *_opti = &opti[0];
//     double *_foptsi = &foptsi[0];
//     double *_comtimei = &comtimei[0];
// 
//     int _n = as<int>(n);
//     int _d = as<int>(d);
//     double _w = as<double>(w);
//     bool _onlyDP = as<bool>(onlyDP);
//     bool _rotated = as<bool>(rotated);
//     bool _isclosed = as<bool>(isclosed);
//     int _skipm = as<int>(skipm);
//     int _autoselectC = as<int>(autoselectC);
//     bool _swap = as<bool>(swap);
// 
// 
//     optimum_reparam(_C1i, _C2i, _n, _d, _w, _onlyDP, _rotated, _isclosed,
//                     _skipm, _autoselectC, _opti, _swap, _foptsi, _comtimei);
// 
//     List ret; ret["opt"] = opti; ret["fopts"] = foptsi;
//     ret["comtime"] = comtimei; ret["swap"] = wrap(_swap);
// 
//     return(ret);
// }
