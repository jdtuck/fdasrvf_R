
#ifndef STIEFEL_H
#define STIEFEL_H

#include "StieVariable.h"
#include "StieVector.h"
#include "Manifold.h"
#include <cmath>
#include "def.h"

// Note that not all the combinations are done.
enum StieMetric{EUCLIDEAN, CANONICAL, STIEMETRICLENGTH};
enum StieRetraction{QF, POLAR, EXP, CONSTRUCTED, STIERETRACTIONLENGTH};
enum StieVectorTransport{PARALLELIZATION, RIGGING, PARALLELTRANSLATION, STIEVECTORTRANSPORTLENGTH};

class Stiefel : public Manifold{
public:
	Stiefel(integer n, integer p);
	virtual ~Stiefel();

	// choose qf, parallelization and intrinsic approach and use householder reflections
	virtual void ChooseStieParamsSet1();

	// choose constructed retraction, parallelization and intrinsic approach
	virtual void ChooseStieParamsSet2();

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

	virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;
	virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

	virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;
	virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

	virtual void CheckParams(void) const;

	virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;
	virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

	virtual void SetParams(PARAMSMAP params);
protected:
	integer n;
	integer p;
	StieMetric metric;
	StieRetraction retraction;
	StieVectorTransport VecTran;

	virtual void ObtainIntrHHR(Variable *x, Vector *etax, Vector *result) const;
	virtual void ObtainExtrHHR(Variable *x, Vector *intretax, Vector *result) const;

	virtual void qfRetraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void qfcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void DiffqfRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

	virtual void ConRetraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void ConcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void DiffConRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

	void ObtainPerp(Variable *x) const;
	virtual void ObtainIntrSquare(Variable *x, Vector *etax, Vector *result) const;
	virtual void ObtainExtrSquare(Variable *x, Vector *intretax, Vector *result) const;
};

#endif
