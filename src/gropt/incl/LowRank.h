#ifndef LOWRANK_H
#define LOWRANK_H

#include <ProductManifold.h>
#include <Stiefel.h>
#include <Euclidean.h>
#include "LowRankVariable.h"
#include "LowRankVector.h"

// Low Rank manifold with representation U \times D times V such that
// X \in \mathcal{M}_r and X = U D V^T, where U \in R^{m \times r}
// D \in R^{r \times r} and V \in R^{n \times r}
class LowRank : public ProductManifold{
public:
	LowRank(integer m, integer n, integer r);
	~LowRank(void);

	//virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

	virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;
	virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

	virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const;
	virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

	virtual void ExtrProjectionStiePerp(Variable *x, Vector *v, Vector *result) const;

	virtual void SetParams(PARAMSMAP params);
protected:
	integer m;
	integer n;
	integer r;
};


#endif // end of LOWRANK_H
