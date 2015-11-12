
#include "LRTRSR1.h"

LRTRSR1::LRTRSR1(const Problem *prob, const Variable *initialx)
{
	const Vector *EMPTYETA;
	if (prob->GetDomain()->GetIsIntrinsic())
		EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
	else
		EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
	Initialization(prob, initialx, EMPTYETA);
	s = EMPTYETA->ConstructEmpty();
	y = EMPTYETA->ConstructEmpty();

	theta = 0.1;
	kappa = 0.9;
	isconvex = false;
	LengthSY = 4;
	S = nullptr;
	Y = nullptr;
	YMGS = nullptr;
	inpss = 0;
	inpsy = 0;
	inpyy = 0;
	Currentlength = 0;
	beginidx = 0;
	SS = nullptr;
	SY = nullptr;
	PMGQ = nullptr;
	P = nullptr;
	gamma = 1;
	ischangedSandY = true;
	SolverName.assign("LRTRSR1");
	prob->SetUseGrad(true);
	prob->SetUseHess(false);
};

LRTRSR1::~LRTRSR1(void)
{
	delete s;
	delete y;
	DeleteVectors(S, LengthSY);
	DeleteVectors(Y, LengthSY);
	DeleteVectors(YMGS, LengthSY);
	if (SS != nullptr)
		delete[] SS;
	if (SY != nullptr)
		delete[] SY;
	if (PMGQ != nullptr)
		delete[] PMGQ;
	if (P != nullptr)
		delete[] P;
};

void LRTRSR1::Run(void)
{
	DeleteVectors(S, LengthSY);
	NewVectors(S, LengthSY);
	DeleteVectors(Y, LengthSY);
	NewVectors(Y, LengthSY);
	DeleteVectors(YMGS, LengthSY);
	NewVectors(YMGS, LengthSY);
	if (SS != nullptr)
		delete[] SS;
	SS = new double[LengthSY * LengthSY];
	if (SY != nullptr)
		delete[] SY;
	SY = new double[LengthSY * LengthSY];
	if (PMGQ != nullptr)
		delete[] PMGQ;
	PMGQ = new double[LengthSY * LengthSY];
	if (P != nullptr)
		delete[] P;
	P = new integer [LengthSY];
	SolversTR::Run();
};

void LRTRSR1::NewVectors(Vector ** &Vs, integer l)
{
	Vs = new Vector *[l];
	for (integer i = 0; i < l; i++)
		Vs[i] = eta1->ConstructEmpty();
};

void LRTRSR1::DeleteVectors(Vector ** &Vs, integer l)
{
	if (Vs != nullptr)
	{
		for (integer i = 0; i < l; i++)
			delete Vs[i];
		delete[] Vs;
	}
};

void LRTRSR1::CheckParams(void)
{
	SolversTR::CheckParams();
	char YES[] = "YES";
	char NO[] = "NO";
	char *status;

	std::cout << "LRTRSR1 METHOD PARAMETERS:" << std::endl;
	status = YES;
	std::cout << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "],\t";
	status = (LengthSY >= 0) ? YES : NO;
	std::cout << "LengthSY      :" << std::setw(15) << LengthSY << "[" << status << "]" << std::endl;
};

void LRTRSR1::HessianEta(Vector *Eta, Vector *result)
{
	integer idx;
	double *v = new double[Currentlength];
	if (ischangedSandY)
	{
		for (integer i = 0; i < Currentlength; i++)
		{
			idx = (i + beginidx) % LengthSY;
			Mani->ScalerVectorAddVector(x1, -gamma, S[idx], Y[idx], YMGS[i]);
		}
		for (integer i = 0; i < Currentlength; i++)
		{
			for (integer j = 0; j < Currentlength; j++)
			{
				PMGQ[i + j * Currentlength] = SY[i + j * LengthSY] - gamma * SS[i + j * LengthSY];
			}
		}
		if (Currentlength > 0)
		{
			// compute LU
			integer info, CurLen = Currentlength;
			dgetrf_(&CurLen, &CurLen, PMGQ, &CurLen, P, &info);
			ischangedSandY = false;
		}
	}

	for (integer i = 0; i < Currentlength; i++)
		v[i] = Mani->Metric(x1, YMGS[i], Eta);

	if (Currentlength > 0)
	{
		char *trans = const_cast<char *> ("n");
		integer info, one = 1, CurLen = Currentlength;
		dgetrs_(trans, &CurLen, &one, PMGQ, &CurLen, P, v, &CurLen, &info);
	}

	Mani->ScaleTimesVector(x1, gamma, Eta, result);
	for (integer i = 0; i < Currentlength; i++)
	{
		Mani->ScalerVectorAddVector(x1, v[i], YMGS[i], result, result);
	}

	delete[] v;
};

void LRTRSR1::UpdateData(void)
{
	double denorminator, norm2ymBs;
	double mintolsq = std::numeric_limits<double>::epsilon();
	double mintol = sqrt(mintolsq);
	Prob->Grad(x2, gf2); ng++;
	eta2->CopyTo(s);
	Mani->InverseVectorTransport(x1, eta2, x2, gf2, eta1); nV++;
	Mani->VectorMinusVector(x1, eta1, gf1, y);
	Mani->VectorMinusVector(x1, y, zeta, zeta);
	denorminator = Mani->Metric(x1, s, zeta);
	inpss = Mani->Metric(x1, s, s);
	norm2ymBs = Mani->Metric(x1, zeta, zeta);
	if (iter == 0) // This is for the robustness when the cost function is quadratic 
	{			   // and its Hessian is identity everywhere.
		inpsy = Mani->Metric(x1, s, y);
		inpyy = Mani->Metric(x1, y, y);
		gamma = inpyy / inpsy;
	}
	if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && norm2ymBs >= mintolsq
		&& (iter != 0 || fabs(gamma - inpsy / inpss) > mintol)) // This is for the robustness when the cost 
												// function is quadratic and its Hessian is identity everywhere.
	{
		inpsy = Mani->Metric(x1, s, y);
		inpyy = Mani->Metric(x1, y, y);
		gamma = inpyy / inpsy;
		if (Currentlength < LengthSY)
		{
			s->CopyTo(S[Currentlength]);
			y->CopyTo(Y[Currentlength]);
			SS[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], S[Currentlength]);
			SY[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[Currentlength]);
			for (integer i = 0; i < Currentlength; i++)
			{
				SS[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], S[i]);
				SS[i + Currentlength * LengthSY] = SS[Currentlength + i * LengthSY];
				SY[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[i]);
				SY[i + Currentlength * LengthSY] = SY[Currentlength + i * LengthSY];
			}
			Currentlength++;
		}
		else
		{
			s->CopyTo(S[beginidx]);
			y->CopyTo(Y[beginidx]);
			for (integer i = 0; i < LengthSY - 1; i++)
			{
				for (integer j = 0; j < LengthSY - 1; j++)
				{
					SS[i + j * LengthSY] = SS[i + 1 + (j + 1) * LengthSY];
					SY[i + j * LengthSY] = SY[i + 1 + (j + 1) * LengthSY];
				}
			}
			SS[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], S[beginidx]);
			SY[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], Y[beginidx]);
			integer idx = 0;
			for (integer i = 0; i < LengthSY - 1; i++)
			{
				idx = (i + beginidx + 1) % LengthSY;
				SS[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, S[idx], S[beginidx]);
				SS[LengthSY - 1 + i * LengthSY] = SS[i + (LengthSY - 1) * LengthSY];
				SY[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, Y[idx], S[beginidx]);
				SY[LengthSY - 1 + i * LengthSY] = SY[i + (LengthSY - 1) * LengthSY];
			}
			beginidx = (++beginidx) % LengthSY;
		}
		isupdated = true;
		ischangedSandY = true;
	}
	else
	{
		isupdated = false;
	}
};

void LRTRSR1::Acceptence(void)
{
	for (integer i = 0; i < Currentlength; i++)
	{
		Mani->VectorTransport(x1, eta2, x2, S[i], S[i]);
		Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]);
	}
	ischangedSandY = true;
};

void LRTRSR1::PrintInfo(void)
{
	printf("\n\tgamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", gamma, inpss, inpsy, inpyy, isupdated);
	std::cout << std::endl;
};

void LRTRSR1::SetParams(PARAMSMAP params)
{
	SolversTR::SetParams(params);
	PARAMSMAP::iterator iter;
	for (iter = params.begin(); iter != params.end(); iter++)
	{
		if (iter->first == static_cast<std::string> ("isconvex"))
		{
			isconvex = ((static_cast<integer> (iter->second)) != 0);
		}
		else
		if (iter->first == static_cast<std::string> ("LengthSY"))
		{
			LengthSY = static_cast<integer> (iter->second);
		}
	}
};
