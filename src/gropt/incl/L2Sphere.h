
#ifndef L2SPHERE_H
#define L2SPHERE_H

#include <Manifold.h>
#include "L2SphereVariable.h"
#include "L2SphereVector.h"
#include <def.h>

enum Repa4NSMetric{ TRAPEZOID, L2SPHEREMETRICLENGTH };
enum Repa4NSRetraction{ NORMALIZED, L2SPHERERETRACTIONLENGTH };
enum Repa4NSVectorTransport{ ISOBYHHR, L2SPHEREVECTORTRANSPORTLENGTH };

// Consider the manifold whose entries are function $l$ satisfying $\int_0^1 l(t)^2 d t = 1$.
// Discretized representation is used.

class L2Sphere : public Manifold{
public:
	L2Sphere(integer inn);
	virtual ~L2Sphere();

	virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;
	virtual void Projection(Variable *x, Vector *v, Vector *result) const;
	virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;
	virtual double Beta(Variable *x, Vector *etax) const;
	virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;
	virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;
	virtual void ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const;

	virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;
	virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

	virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;
	virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

	virtual void CheckParams(void) const;

	virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;
	virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;
private:
	mutable integer n;

	Repa4NSMetric metric;
	Repa4NSRetraction retraction;
	Repa4NSVectorTransport VecTran;
};

#endif // end of L2SPHERE_H
