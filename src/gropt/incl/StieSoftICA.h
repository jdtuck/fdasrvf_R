
#ifndef STIESOFTICA_H
#define STIESOFTICA_H

#include "Stiefel.h"
#include "StieVariable.h"
#include "StieVector.h"
#include "Problem.h"
#include "SharedSpace.h"
#include "def.h"

// (12.1.1) in Wen Huang's Thesis
class StieSoftICA : public Problem{
public:
	StieSoftICA(double *inCs, integer inn, integer inp, integer inN);
	virtual double f(Variable *x) const;

	virtual void EucGrad(Variable *x, Vector *egf) const;
	virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	double *Cs;
	mutable integer n;
	mutable integer p;
	mutable integer N;
};

#endif // end of STIESOFTICA_H
