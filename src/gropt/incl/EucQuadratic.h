
#ifndef EUCQUADRATIC_H
#define EUCQUADRATIC_H

#include "Euclidean.h"
#include "EucVariable.h"
#include "EucVector.h"
#include "Problem.h"
#include "SharedSpace.h"
#include "def.h"

// min_x x^T A x, where A is a symmetric positive definite matrix
class EucQuadratic : public Problem{
public:
	EucQuadratic(double *M, integer dim);
	virtual ~EucQuadratic();
	virtual double f(Variable *x) const;
	virtual void Grad(Variable *x, Vector *gf) const;
	virtual void HessianEta(Variable *x, Vector *etax, Vector *xix) const;

	double *A;
	integer Dim;
};

#endif // end of EUCQUADRATIC_H
