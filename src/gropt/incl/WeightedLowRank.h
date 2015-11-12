
#ifndef WEIGHTEDLOWRANK_H
#define WEIGHTEDLOWRANK_H

#include "Stiefel.h"
#include "StieVariable.h"
#include "StieVector.h"
#include "Problem.h"
#include "SharedSpace.h"
#include "def.h"
#include "Element.h"
#include "ProductElement.h"
#include "ProductManifold.h"

class WeightedLowRank : public Problem{
public:
	WeightedLowRank(double *inA, double *inW, integer inm, integer inn, integer inr);
	virtual ~WeightedLowRank();
	virtual double f(Variable *x) const;

	//virtual void EucGrad(Variable *x, Vector *egf) const;
	//virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	virtual void RieGrad(Variable *x, Vector *gf) const;
	virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;
	
	double *A;
	double *W;
	integer m;
	integer n;
	integer r;

};
#endif
