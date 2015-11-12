
#ifndef STIESUMBROCKETT_H
#define STIESUMBROCKETT_H

#include "Stiefel.h"
#include "StieVariable.h"
#include "StieVector.h"
#include <ProductElement.h>
#include <ProductManifold.h>
#include "Problem.h"
#include "SharedSpace.h"
#include "def.h"

// min_X X^T B X D, where B is a symmetric positive definite matrix, D is a diagonal matrix
// and X \in St(p, n).
class StieSumBrockett : public Problem{
public:
	StieSumBrockett(double *inB1, double *inD1, double *inB2, double *inD2, double *inB3, 
				 double *inD3, integer inn, integer inp, integer inm, integer inq);
	virtual ~StieSumBrockett();
	virtual double f(Variable *x) const;

	virtual void EucGrad(Variable *x, Vector *egf) const;
	virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	double *B1;
	double *D1;
	double *B2;
	double *D2;
	double *B3;
	double *D3;
	integer n;
	integer p;
	integer m;
	integer q;
};

#endif // end of STIESUMBROCKETT_H
