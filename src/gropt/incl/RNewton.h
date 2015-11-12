#ifndef RNEWTON_H
#define RNEWTON_H

#include "SolversLS.h"
#include "def.h"

enum tCGLSstatusSet{ LS_NEGCURVTURE, LS_LCON, LS_SCON, LS_MAXITER, TCGLSSTATUSSETLENGTH };

class RNewton : public SolversLS{
public:
	RNewton(const Problem *prob, const Variable *initialx);
	virtual void CheckParams();
	virtual ~RNewton();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool useRand;
	integer Max_Inner_Iter;
	integer Min_Inner_Iter;
	double theta;
	double kappa;
protected:
	virtual void GetSearchDir();
	void PreConditioner(Variable *x, Vector *eta, Vector *result);
	void HessianEta(Vector *Eta, Vector *result) const;
	virtual void PrintInfo();
	void tCG_LS();

	Vector *r, *z, *delta, *Hd;
	integer innerIter;
	tCGLSstatusSet tCGLSstatus;
	std::string *tCGLSstatusSetnames;
};

#endif // end of RNEWTON_H
