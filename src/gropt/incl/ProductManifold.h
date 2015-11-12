#ifndef PRODUCTMANIFOLD_H
#define PRODUCTMANIFOLD_H

#define CHECKMEMORY

#define ProdVariable ProductElement
#define ProdVector ProductElement

#include "Manifold.h"
#include "ProductElement.h"
#include <EucVariable.h>
#include <EucVector.h>
#include <cstdarg>
#include <map>
#include "def.h"

class ProductManifold : public Manifold{
public:
	ProductManifold(integer numberofmanifolds, ...);
	ProductManifold(Manifold **inmanifolds, integer innumofmani, integer *inpowsinterval, integer innumoftotalmani);
	virtual ~ProductManifold();

	virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;
	virtual void LinearOPEEta(Variable *x, LinearOPE *Hx, Vector *etax, Vector *result) const;
	virtual void ScaleTimesVector(Variable *x, double scaler, Vector *etax, Vector *result) const;
	virtual void VectorLinearCombination(Variable *x, double scaler1, Vector *etax, double scaler2, Vector *xix, Vector *result) const;
	virtual void Projection(Variable *x, Vector *v, Vector *result) const;
	virtual void RandomTangentVectors(Variable *x, integer N, Vector **result_arr) const;
	virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;
	virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;
	virtual double Beta(Variable *x, Vector *etax) const;
	virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;
	virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;
	virtual void HaddScaledRank1OPE(Variable *x, LinearOPE *Hx, double scaler, Vector *etax, Vector *xix, LinearOPE *result) const;

	virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

	virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;
	virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

	virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;
	virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

	virtual void CheckParams(void) const;

	inline Manifold *GetManifold(integer idx) const { if (idx < numofmani) return manifolds[idx]; return nullptr; };

protected:
	void ProdElementToElement(const ProductElement *ProdElem, Element *Elem) const;
	void ElementToProdElement(const Element *Elem, ProductElement *ProdElem) const;

	Manifold **manifolds;
	integer numofmani;
	integer *powsinterval;
	integer numoftotalmani;
};

#endif // end of PRODUCTMANIFOLD_H