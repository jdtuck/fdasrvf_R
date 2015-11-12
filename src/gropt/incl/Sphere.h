#ifndef SPHERE_H
#define SPHERE_H

#include "SphereVariable.h"
#include "SphereVector.h"
#include <Stiefel.h>

class Sphere : public Stiefel{
public:
	Sphere(integer n);
	virtual ~Sphere();

	// choose exponential map, parallel translation and extrinsic approach and no householder reflections
	// Even though the Householder reflections are not used, the locking condition is satisfied.
	void ChooseSphereParamsSet1();

	// choose qf, parallel translation and extrinsic approach and no householder reflections
	// The locking conidition is not satisfied
	void ChooseSphereParamsSet2();

	// choose qf, parallel translation and extrinsic approach and no householder reflections
	// Beta \neq 1 is used and the locking conidition is satisfied
	void ChooseSphereParamsSet3();

	virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;
	virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;
	virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

	virtual void SetParams(PARAMSMAP params);
private:
	// functions for exponential map and parallel translation using extrinsic approach
	virtual void ExpRetraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void ExpcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void ExpDiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;
	virtual void ExpVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;
	virtual void ExpInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void ExpHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void ExpTranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void ExpTranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;
};

#endif // end of SPHERE_H