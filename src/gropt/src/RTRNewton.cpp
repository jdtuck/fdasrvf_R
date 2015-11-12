
#include "RTRNewton.h"

RTRNewton::RTRNewton(const Problem *prob, const Variable *initialx)
{
	const Vector *EMPTYETA;
	if (prob->GetDomain()->GetIsIntrinsic())
		EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
	else
		EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
	SolversTR::Initialization(prob, initialx, EMPTYETA);
	SolverName.assign("RTRNewton");
	prob->SetUseGrad(true);
	prob->SetUseHess(true);
};

void RTRNewton::HessianEta(Vector *Eta, Vector *result)
{
	Prob->HessianEta(x1, Eta, result);
};
