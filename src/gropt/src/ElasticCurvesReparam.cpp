
#include "ElasticCurvesReparam.h"
#include <Rcpp.h>
using namespace Rcpp;

std::map<integer *, integer> *CheckMemoryDeleted;

RcppExport SEXP opt_reparam(SEXP C1, SEXP C2, SEXP n, SEXP d, SEXP w,
                            SEXP onlyDP, SEXP rotated, SEXP isclosed,
                            SEXP skipm, SEXP autoselectC, SEXP opt, SEXP swap,
                            SEXP fopts, SEXP comtime){

    NumericVector C1i(C1);
    NumericVector C2i(C2);
    NumericVector opti(opt);
    NumericVector foptsi(fopts);
    NumericVector comtimei(comtime);

    double *_C1i = &C1i[0];
    double *_C2i = &C2i[0];
    double *_opti = &opti[0];
    double *_foptsi = &foptsi[0];
    double *_comtimei = &comtimei[0];

    int _n = as<int>(n);
    int _d = as<int>(d);
    double _w = as<double>(w);
    bool _onlyDP = as<bool>(onlyDP);
    bool _rotated = as<bool>(rotated);
    bool _isclosed = as<bool>(isclosed);
    int _skipm = as<int>(skipm);
    int _autoselectC = as<int>(autoselectC);
    bool _swap = as<bool>(swap);


    optimum_reparam(_C1i, _C2i, _n, _d, _w, _onlyDP, _rotated, _isclosed,
                    _skipm, _autoselectC, _opti, _swap, _foptsi, _comtimei);

    List ret; ret["opt"] = opti; ret["fopts"] = foptsi;
    ret["comtime"] = comtimei; ret["swap"] = wrap(_swap);

    return(ret);
}

void optimum_reparam(double *C1, double *C2, int n, int d, double w,
        bool onlyDP, bool rotated, bool isclosed, int skipm, int autoselectC,
        double *opt, bool swap, double *fopts, double *comtime)
{
    /* dimensions of input matrices */
    /* opt size is n + d*d +1 */
    /* fopts and comtime are 5 x 1*/
    integer n1, d1;
    n1 = static_cast<integer> (n);
    d1 = static_cast<integer> (d);
    bool swapi;

    std::string methodname = "";
    if (!onlyDP)
        methodname = "LRBFGS";

    init_genrand(0);

    CheckMemoryDeleted = new std::map<integer *, integer>;

    integer numofmanis = 3;
    integer numofmani1 = 1;
    integer numofmani2 = 1;
    integer numofmani3 = 1;
    L2SphereVariable FNSV(n);
    OrthGroupVariable OGV(d);
    EucVariable EucV(1);
    ProductElement *Xopt = new ProductElement(numofmanis, &FNSV, numofmani1, &OGV, numofmani2, &EucV, numofmani3);

    integer ns, lms;

    DriverElasticCurvesRO(C1, C2, d1, n1, w, rotated, isclosed, onlyDP, skipm, methodname,
            autoselectC, Xopt, swapi, fopts, comtime, ns, lms);

    swap = swapi;

    /* get output data */
    integer sizex = n1 + d1 * d1 + 1;
    const double *Xoptptr = Xopt->ObtainReadData();
    integer inc = 1;
    dcopy_(&sizex, const_cast<double *> (Xoptptr), &inc, opt, &inc);

    delete Xopt;

    std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
    for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
    {
        if (iter->second != 1)
            std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
    }
    delete CheckMemoryDeleted;
    return;
}
