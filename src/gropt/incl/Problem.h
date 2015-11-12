#ifndef PROBLEM_H
#define PROBLEM_H

#include <Element.h>
#include "Manifold.h"
#include <cmath>
#include <iostream>
#include "def.h"

class Manifold;

class Problem{
public:
	virtual ~Problem(void) = 0;
	virtual double f(Variable *x) const = 0;
	virtual void Grad(Variable *x, Vector *gf) const;
	virtual void HessianEta(Variable *x, Vector *etax, Vector *xix) const;

	virtual void RieGrad(Variable *x, Vector *gf) const;

	// etax and xix must be different arguments
	virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

	virtual void EucGrad(Variable *x, Vector *egf) const;

	// etax and xix must be different arguments
	virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	virtual void CheckGradHessian(const Variable *x) const;

	virtual void SetDomain(Manifold *inDomain);
	inline Manifold *GetDomain(void) const { return Domain; };

	virtual void SetUseGrad(bool usegrad) const;
	virtual void SetUseHess(bool usehess) const;
	inline bool GetUseGrad(void) const { return UseGrad; };
	inline bool GetUseHess(void) const { return UseHess; };
protected:
	Manifold *Domain;
	mutable bool UseGrad;
	mutable bool UseHess;
};

#endif // end of PROBLEM_H
