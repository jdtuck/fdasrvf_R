#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H

#include "ForDebug.h"
#include "EucVariable.h"
#include "EucVector.h"
#include "Manifold.h"
#include "def.h"

class Euclidean : public Manifold{
public:
	Euclidean(integer r, integer c = 1, integer n = 1);
	virtual ~Euclidean();
	virtual void CheckParams(void) const;

	virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;
	virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

	integer row;
	integer col;
	integer num;
};

#endif // end of EUCLIDEAN_H
