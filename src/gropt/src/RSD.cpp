
#include "RSD.h"

RSD::RSD(const Problem *prob, const Variable *initialx)
{
	const Vector *EMPTYETA;
	if (prob->GetDomain()->GetIsIntrinsic())
		EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
	else
		EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
	Initialization(prob, initialx, EMPTYETA);
	SolverName.assign("RSD");
	prob->SetUseGrad(true);
	prob->SetUseHess(false);
};

void RSD::GetSearchDir(void)
{
	Mani->ScaleTimesVector(x1, -1, gf1, eta1);
};
