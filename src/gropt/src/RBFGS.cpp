
#include "RBFGS.h"

RBFGS::RBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
{
	const Vector *EMPTYETA;
	if (prob->GetDomain()->GetIsIntrinsic())
		EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
	else
		EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
	Initialization(prob, initialx, EMPTYETA);
	bool initHisnull = (initialH == nullptr);
	if (initHisnull)
	{
		if (prob->GetDomain()->GetIsIntrinsic())
		{
			initialH = new LinearOPE(prob->GetDomain()->GetEMPTYINTR()->Getlength());
		}
		else
		{
			initialH = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR()->Getlength());
		}
	}
	initialH->ScaledIdOPE();

	H = initialH->ConstructEmpty();
	tildeH = initialH->ConstructEmpty();
	initialH->CopyTo(H);
	s = EMPTYETA->ConstructEmpty();
	y = EMPTYETA->ConstructEmpty();

	isconvex = false;
	nu = 1e-4;
	mu = 1;
	SolverName.assign("RBFGS");
	if (initHisnull)
		delete initialH;
	prob->SetUseGrad(true);
	prob->SetUseHess(false);
};

RBFGS::~RBFGS(void)
{
	delete s;
	delete y;
	delete H;
	delete tildeH;
};

void RBFGS::CheckParams(void)
{
	SolversLS::CheckParams();
	char YES[] = "YES";
	char NO[] = "NO";
	char *status;

	std::cout << "RBFGS METHOD PARAMETERS:" << std::endl;
	status = (nu >= 0 && nu < 1) ? YES : NO;
	std::cout << "nu            :" << std::setw(15) << nu << "[" << status << "],\t";
	status = (mu >= 0) ? YES : NO;
	std::cout << "mu            :" << std::setw(15) << mu << "[" << status << "]" << std::endl;
	status = YES;
	std::cout << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "]" << std::endl;
};

void RBFGS::GetSearchDir(void)
{
	Mani->LinearOPEEta(x1, H, gf1, eta1); nH++;
	Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
};

void RBFGS::UpdateData(void)
{
	Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
	Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
	betay = Mani->Beta(x1, eta2);
	Mani->ScalerVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
	inpsy = Mani->Metric(x2, s, y);
	if (isconvex && iter == 1 && inpsy > 0)
		H->ScaledIdOPE(inpsy / Mani->Metric(x2, y, y));

	Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
	inpss = Mani->Metric(x2, s, s);
	if (inpsy / inpss >= nu * pow(ngf, mu) && inpss > std::numeric_limits<double>::epsilon() 
        && inpsy > std::numeric_limits<double>::epsilon())
	{
		Mani->LinearOPEEta(x2, tildeH, y, zeta); // zeta = tildeH y
		Mani->HaddScaledRank1OPE(x2, tildeH, -1.0 / inpsy, s, zeta, H);
		Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y

		Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, zeta, s, H);
		Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
		isupdated = true;
	}
	else
	{
		isupdated = false;
		tildeH->CopyTo(H);
#ifdef TESTELASTICCURVESRO
		H->ScaledIdOPE(1);
#endif
	}
};

void RBFGS::PrintInfo(void)
{
	printf("\n\tbetay:%.3e,inpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", betay, inpss, inpsy, isupdated);
	std::cout << std::endl;
};

void RBFGS::SetParams(PARAMSMAP params)
{
	SolversLS::SetParams(params);
	PARAMSMAP::iterator iter;
	for (iter = params.begin(); iter != params.end(); iter++)
	{
		if (iter->first == static_cast<std::string> ("isconvex"))
		{
			isconvex = ((static_cast<integer> (iter->second)) != 0);
		}
		else
		if (iter->first == static_cast<std::string> ("nu"))
		{
			nu = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("mu"))
		{
			mu = iter->second;
		}
	}
};
