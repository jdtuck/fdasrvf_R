#ifndef MANIFOLD_H
#define MANIFOLD_H

#define Variable Element
#define Vector Element

#include <Problem.h>
#include "Element.h"
#include "LinearOPE.h"
#include "SharedSpace.h"
#include <iomanip>
#include <cmath>
#include "def.h"

class Problem;
class ProductManifold;

class Manifold{
public:
    virtual ~Manifold(void) = 0;
	virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

	// etax and result must be different arguments
	virtual void LinearOPEEta(Variable *x, LinearOPE *Hx, Vector *etax, Vector *result) const;
	
	// etax and result can be a same argument
	virtual void ScaleTimesVector(Variable *x, double scaler, Vector *etax, Vector *result) const;

	// (etax and result) or (xix and result) can be a same argument
	virtual void VectorAddVector(Variable *x, Vector *etax, Vector *xix, Vector *result) const;

	// (etax and result) or (xix and result) can be a same argument
	virtual void VectorMinusVector(Variable *x, Vector *etax, Vector *xix, Vector *result) const;

	// (etax and result) or (xix and result) can be a same argument
	virtual void ScalerVectorAddVector(Variable *x, double scaler, Vector *etax, Vector *xix, Vector *result) const;

	// (etax and result) or (xix and result) can be a same argument
	virtual void ScalerVectorMinusVector(Variable *x, double scaler, Vector *etax, Vector *xix, Vector *result) const;

	// (etax and result) or (xix and result) can be a same argument
	virtual void VectorLinearCombination(Variable *x, double scaler1, Vector *etax, double scaler2, Vector *xix, Vector *result) const;

	// etax and result can be a same argument
	virtual void Projection(Variable *x, Vector *etax, Vector *result) const;

	// for Gradient sampling method. not used now.
	virtual void RandomTangentVectors(Variable *x, integer N, Vector **result_arr) const; 

	virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

	// xiy and result can be a same argument
	virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

	// xix and result can be a same argument
	virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

	virtual double Beta(Variable *x, Vector *etax) const;

	// xix and result can be a same argument
	virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

	// xiy and result can be a same argument
	virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

	// Hx and result can be a same argument
	virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

	// Hx and result can be a same argument
	virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

	// Hx and result can be a same argument
	virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

	// Hx and result can be a same argument
	virtual void HaddScaledRank1OPE(Variable *x, LinearOPE *Hx, double scaler, Vector *etax, Vector *xix, LinearOPE *result) const;

	// etax and etaxflat can be a same argument
	virtual void ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const;

	// etax and result must be different arguments
	virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;
	// intretax and result must be different arguments
	virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

	// etax and result can be a same argument
	virtual void IntrProjection(Variable *x, Vector *etax, Vector *result) const;

	// etax and result can be a same argument
	virtual void ExtrProjection(Variable *x, Vector *etax, Vector *result) const;

	inline std::string GetName(void) const { return name; };
	inline integer GetIntrDim(void) const { return IntrinsicDim; };
	inline integer GetExtrDim(void) const { return ExtrinsicDim; };
	inline bool GetIsIntrinsic(void) const { return IsIntrApproach; };
	inline const Vector *GetEMPTYINTR(void) const { return EMPTYINTR; };
	inline const Vector *GetEMPTYEXTR(void) const { return EMPTYEXTR; };
	inline bool GetHasLockCon(void) const { return HasLockCon; };

	virtual void CheckParams(void) const;

	// egf and result can be a same argument
	virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *result, const Problem *prob) const = 0;

	// exix and result can be a same argument
	virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const = 0;

	virtual void SetParams(PARAMSMAP params);

	virtual void CheckIntrExtr(Variable *x) const;
	virtual void CheckRetraction(Variable *x) const;
	virtual void CheckDiffRetraction(Variable *x, bool IsEtaXiSameDir = true) const;
	virtual void CheckLockingCondition(Variable *x) const;
	virtual void CheckcoTangentVector(Variable *x) const;
	virtual void CheckIsometryofVectorTransport(Variable *x) const;
	virtual void CheckIsometryofInvVectorTransport(Variable *x) const;
	virtual void CheckVecTranComposeInverseVecTran(Variable *x) const;
	virtual void CheckTranHInvTran(Variable *x) const;
	virtual void CheckHaddScaledRank1OPE(Variable *x) const;

	inline void SetHasHHR(bool inHasHHR) { HasHHR = inHasHHR; };
	inline bool GetHasHHR(void) { return HasHHR; };
	inline void SetIsIntrApproach(bool IsIntrAppr){ IsIntrApproach = IsIntrAppr; };
protected:
	mutable bool HasHHR;
	bool UpdBetaAlone;
	std::string name;
	bool IsIntrApproach;
	// metric
	integer IntrinsicDim;
	integer ExtrinsicDim;
	Vector *EMPTYINTR;
	Vector *EMPTYEXTR;

	bool HasLockCon;

	// The next 5 functions modified above 5 functions such that the locking condition is satisfied.
	// They are not required to be satisfied.
	virtual void LCVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;
	virtual void LCInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;
	virtual void LCHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void LCTranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;
	virtual void LCTranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

	virtual void Obtainnu1nu2forLC(Variable *x, Vector *etax, Variable *y) const;
};

#endif // end of MANIFOLD_H