#ifndef MEXPROBLEM_H
#define MEXPROBLEM_H

#include "Problem.h"
#include <cstring>
#include <def.h>

#ifdef MATLAB_MEX_FILE

class mexProblem : public Problem{
public:
	mexProblem(const mxArray *inf, const mxArray *ingf, const mxArray *inHess);
	virtual ~mexProblem();
	virtual double f(Variable *x) const;

	virtual void EucGrad(Variable *x, Vector *egf) const;
	virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	static void ObtainMxArrayFromElement(mxArray *&Xmx, const Element *X);
	static void ObtainElementFromMxArray(Element *X, const mxArray *Xmx);
	static mxArray *GetFieldbyName(const mxArray *S, integer idxstruct, const char *name);
protected:
	const mxArray *mxf;
	const mxArray *mxgf;
	const mxArray *mxHess;
};

#endif // end of MATLAB_MEX_FILE

#endif // end of MEXPROBLEM_H
