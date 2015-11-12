
#include "RTRSD.h"

RTRSD::RTRSD(const Problem *prob, const Variable *initialx)
{
	const Vector *EMPTYETA;
	if (prob->GetDomain()->GetIsIntrinsic())
		EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
	else
		EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
	SolversTR::Initialization(prob, initialx, EMPTYETA);
	theta = 0.1;
	kappa = 0.9;
	SolverName.assign("RTRSD");
	prob->SetUseGrad(true);
	prob->SetUseHess(false);
};

void RTRSD::HessianEta(Vector *Eta, Vector *result)
{
	Eta->CopyTo(result);
};
