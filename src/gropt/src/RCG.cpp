
#include "RCG.h"

RCG::RCG(const Problem *prob, const Variable *initialx)
{
	const Vector *EMPTYETA;
	if (prob->GetDomain()->GetIsIntrinsic())
		EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
	else
		EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
	Initialization(prob, initialx, EMPTYETA);
	sigma = 0;
	RCGmethod = POLAK_RIBIERE_MOD;
	ManDim = std::numeric_limits<integer>::max();

	SolverName.assign("RCG");
	prob->SetUseGrad(true);
	prob->SetUseHess(false);

	RCGmethodSetnames = new std::string [RCGMETHODSLENGTH];
	RCGmethodSetnames[FLETCHER_REEVES].assign("FLETCHER_REEVES");
	RCGmethodSetnames[POLAK_RIBIERE_MOD].assign("POLAK_RIBIERE_MOD");
	RCGmethodSetnames[HESTENES_STIEFEL].assign("HESTENES_STIEFEL");
	RCGmethodSetnames[FR_PR].assign("FR_PR");
	RCGmethodSetnames[DAI_YUAN].assign("DAI_YUAN");
	RCGmethodSetnames[HAGER_ZHANG].assign("HAGER_ZHANG");
};

RCG::~RCG(void)
{
	delete[] RCGmethodSetnames;
};

void RCG::CheckParams(void)
{
	SolversLS::CheckParams();
	char YES[] = "YES";
	char NO[] = "NO";
	char *status;

	std::cout << "RCG METHOD PARAMETERS:" << std::endl;
	status = (ManDim >= 0 && ManDim <= std::numeric_limits<integer>::max()) ? YES : NO;
	std::cout << "ManDim        :" << std::setw(15) << ManDim << "[" << status << "],\t";
	status = (RCGmethod >= 0 && RCGmethod <= RCGMETHODSLENGTH) ? YES : NO;
	std::cout << "RCGmethod     :" << std::setw(15) << RCGmethodSetnames[RCGmethod] << "[" << status << "]" << std::endl;
};

void RCG::PrintInfo(void)
{
	if (iter % ManDim == 0 || Mani->Metric(x1, eta1, gf1) >= -std::numeric_limits<double>::epsilon()) // restart and safeguard
		printf("\n\tsigma:%.3e,Reset search direction to the negative gradient,", sigma);
	else
		printf("\n\tsigma:%.3e,", sigma);
	std::cout << std::endl;
};

void RCG::GetSearchDir(void)
{
	if (iter % ManDim == 0 || Mani->Metric(x1, eta1, gf1) >= -std::numeric_limits<double>::epsilon()) // restart and safeguard
	{
		Mani->ScaleTimesVector(x1, -1.0, gf1, eta1);
	}
};

void RCG::UpdateData(void)
{
	if (iter % ManDim != 0)
	{
		if (RCGmethod == FLETCHER_REEVES)
		{
			sigma = Mani->Metric(x2, gf2, gf2) / Mani->Metric(x1, gf1, gf1);
			Mani->VectorTransport(x1, eta2, x2, eta1, zeta);
		}
		else
		if (RCGmethod == POLAK_RIBIERE_MOD)
		{
			Mani->VectorTransport(x1, eta2, x2, gf1, zeta);
			Mani->VectorMinusVector(x2, gf2, zeta, zeta);
			sigma = Mani->Metric(x2, zeta, gf2) / Mani->Metric(x1, gf1, gf1);
			if (LineSearch_LS == STRONGWOLFE && sigma <= 0)
				sigma = 0;
			else
				Mani->VectorTransport(x1, eta2, x2, eta1, zeta);
		}
		else
		if (RCGmethod == HESTENES_STIEFEL)
		{
			double numerator, denominator;
			Mani->VectorTransport(x1, eta2, x2, eta1, zeta);

			Mani->VectorTransport(x1, eta2, x2, gf1, eta1);
			Mani->VectorMinusVector(x2, gf2, eta1, eta1);
			numerator = Mani->Metric(x2, eta1, gf2);
			denominator = Mani->Metric(x2, zeta, eta1);
			sigma = numerator / denominator;
		}
		else
		if (RCGmethod == FR_PR)
		{
			double sigmaFR, sigmaPR;
			Mani->VectorTransport(x1, eta2, x2, gf1, zeta);
			Mani->VectorMinusVector(x2, gf2, zeta, zeta);
			sigmaPR = Mani->Metric(x2, zeta, gf2) / Mani->Metric(x1, gf1, gf1);
			sigmaFR = Mani->Metric(x2, gf2, gf2) / Mani->Metric(x1, gf1, gf1);
			if (sigmaPR < -sigmaFR)
				sigma = -sigmaFR;
			else
			if (sigmaPR > sigmaFR)
				sigma = sigmaFR;
			else
				sigma = sigmaPR;
			Mani->VectorTransport(x1, eta2, x2, eta1, zeta);
		}
		else
		if (RCGmethod == DAI_YUAN)
		{
			Mani->VectorTransport(x1, eta2, x2, eta1, zeta);

			Mani->VectorTransport(x1, eta2, x2, gf1, eta1);
			Mani->VectorMinusVector(x2, gf2, eta1, eta1);
			sigma = Mani->Metric(x2, gf2, gf2) / Mani->Metric(x2, zeta, eta1);
		}
		else
		if (RCGmethod == HAGER_ZHANG)
		{
			double temp1, temp2;
			Mani->VectorTransport(x1, eta2, x2, eta1, zeta);

			Mani->VectorTransport(x1, eta2, x2, gf1, eta1);
			Mani->VectorMinusVector(x2, gf2, eta1, eta1);
			temp1 = Mani->Metric(x2, eta1, zeta);
			temp2 = - 2.0 * Mani->Metric(x2, eta1, eta1) / temp1;
			Mani->ScalerVectorAddVector(x2, temp2, zeta, eta1, eta1);
			sigma = Mani->Metric(x2, eta1, gf2) / temp1;
		}
		Mani->ScalerVectorMinusVector(x2, sigma, zeta, gf2, eta1);
	}
};

void RCG::SetParams(PARAMSMAP params)
{
	SolversLS::SetParams(params);
	PARAMSMAP::iterator iter;
	for (iter = params.begin(); iter != params.end(); iter++)
	{
		if (iter->first == static_cast<std::string> ("ManDim"))
		{
			ManDim = static_cast<integer> (iter->second);
		}
		else
		if (iter->first == static_cast<std::string> ("RCGmethod"))
		{
			RCGmethod = static_cast<RCGmethods> (static_cast<integer> (iter->second));
		}
	}
};
