
#include "SolversTR.h"

void SolversTR::Run(void)
{
	Variable *xTemp;
	Vector *gfTemp;
	starttime = getTickCount();
	Solvers::Run();
	double sqeps = sqrt(std::numeric_limits<double>::epsilon());

	f1 = Prob->f(x1);
	Prob->Grad(x1, gf1);

	ngf0 = sqrt(Mani->Metric(x1, gf1, gf1));
	ngf = ngf0;

	iter = 0;
	if (DEBUG >= ITERRESULT)
	{
		printf("i:%d,f:%.3e,|gf|:%.3e,\n", iter, f1, ngf);
		timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
		funSeries[iter] = f1; gradSeries[iter] = ngf;
		nf++; ng++;
	}
	Delta = initial_Delta;
	bool isstop = false;
	while (((! isstop) && iter < Max_Iteration) || iter < Min_Iteration)
	{
		InitialVector(); // Obtain initial guess, eta1, for local model
		tCG_TR(); // obtain eta2
		Mani->Retraction(x1, eta2, x2);	nR++;
		f2 = Prob->f(x2); nf++;
		HessianEta(eta2, zeta); nH++; // Hessian * eta2

		Mani->ScalerVectorAddVector(x1, 0.5, zeta, gf1, eta1);
		rho = (f1 - f2) / (-Mani->Metric(x1, eta2, eta1));
		UpdateData(); // Update S&Y or H or B

		if (rho > 0.75)
		{
			if (tCGstatus == TR_EXCREGION || tCGstatus == TR_NEGCURVTURE)
				Delta *= Magnified_tau;
			if (Delta > maximum_Delta)
			{
				std::cout << "reach the maximum of radius" << std::endl;
				Delta = maximum_Delta;
			}
		}
		else
		if (rho < 0.25)
		{
			Delta *= Shrinked_tau;
			if (Delta < minimum_Delta)
			{
				std::cout << "reach the minimum of radius" << std::endl;
				break;
			}
		}
		isstop = IsStopped();
		if (rho > Acceptence_Rho || (fabs(f1 - f2) / (fabs(f1) + 1) < sqeps && f2 < f1))
		{
			Acceptence(); // Algorithm specific operations
			ngf = sqrt(Mani->Metric(x2, gf2, gf2));
			xTemp = x1; x1 = x2; x2 = xTemp;
			gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
            iter++;
			if (DEBUG >= ITERRESULT && iter % OutputGap == 0)
			{
				PrintGenInfo();
				PrintInfo(); // Output information specific to Algorithms
			}
			f1 = f2;
		}
		else
		{
            iter++;
			if (DEBUG >= ITERRESULT && iter % OutputGap == 0)
			{
				std::cout << "X_{" << iter << "} WAS REJECTED." << std::endl;
				PrintGenInfo();
				PrintInfo(); // Output information specific to Algorithms
			}
		}

		if (DEBUG >= ITERRESULT)
		{
			timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
			funSeries[iter] = f2; gradSeries[iter] = ngf;
		}
	}
	ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
	if (DEBUG >= ITERRESULT)
		lengthSeries = iter + 1;

	if (DEBUG >= FINALRESULT)
	{
		printf("Iter:%d,f:%.3e,|gf|:%.3e,|gf|/|gf0|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,", iter, f2,
			ngf, ngf / ngf0, ComTime, nf, ng, nR);
		if (nH != 0)
		{
			printf("nH:%d,", nH);
		}
		if (nV != 0)
		{
			printf("nV(nVp):%d(%d),", nV, nVp);
		}
		printf("\n");
	}
};

void SolversTR::InitialVector(void)
{
	Mani->ScaleTimesVector(x1, 0, gf1, eta1);
};

void SolversTR::tCG_TR(void)
{
	double e_Pe, r_r, norm_r, norm_r0, d_Pd, z_r, e_Pd, d_Hd, alphatemp, e_Pe_new, tempnum, zold_rold, betatemp, tautemp;
	integer j;

	if (useRand)
	{
		HessianEta(eta1, r); nH++;
		Mani->VectorAddVector(x1, gf1, r, r);
		e_Pe = Mani->Metric(x1, eta1, eta1);
	}
	else
	{
		gf1->CopyTo(r);
		e_Pe = 0;
	}

	r_r = Mani->Metric(x1, r, r);
	norm_r = sqrt(r_r);
	norm_r0 = norm_r;

	PreConditioner(x1, r, z);
	
	z_r = Mani->Metric(x1, z, r);
	d_Pd = z_r;

	Mani->ScaleTimesVector(x1, -1.0, z, delta);

	if (useRand)
		e_Pd = Mani->Metric(x1, eta1, delta);
	else
		e_Pd = 0;

	tCGstatus = TR_MAXITER; /* pre-assume termination j == max_inner*/

	eta1->CopyTo(eta2);

	for (j = 0; j < Max_Inner_Iter; j++)
	{
		HessianEta(delta, Hd); nH++;
		d_Hd = Mani->Metric(x1, delta, Hd);
		alphatemp = z_r / d_Hd;
		e_Pe_new = e_Pe + 2.0 * alphatemp * e_Pd + alphatemp * alphatemp * d_Pd;

		if (d_Hd <= 0 || e_Pe_new >= (Delta * Delta))
		{
			tautemp = (-e_Pd + sqrt(e_Pd * e_Pd + d_Pd * (Delta * Delta - e_Pe))) / d_Pd;
			Mani->ScalerVectorAddVector(x1, tautemp, delta, eta2, eta2);

			if (d_Hd < 0)
				tCGstatus = TR_NEGCURVTURE; /* negative curvature*/
			else
				tCGstatus = TR_EXCREGION; /* exceeded trust region*/
			break;
		}
		e_Pe = e_Pe_new;
		Mani->ScalerVectorAddVector(x1, alphatemp, delta, eta2, eta2);

		Mani->ScalerVectorAddVector(x1, alphatemp, Hd, r, r);

		Mani->Projection(x1, r, zeta);
		zeta->CopyTo(r);

		r_r = Mani->Metric(x1, r, r);
		norm_r = sqrt(r_r);

		tempnum = pow(norm_r0, theta);
		if (j >= Min_Inner_Iter && norm_r <= norm_r0 * ((tempnum < kappa) ? tempnum : kappa))
		{
			if (kappa < tempnum)
				tCGstatus = TR_LCON; /* linear convergence*/
			else
				tCGstatus = TR_SCON; /* superlinear convergence*/
			break;
		}

		PreConditioner(x1, r, z);

		zold_rold = z_r;
		z_r = Mani->Metric(x1, z, r);

		betatemp = z_r / zold_rold;
		Mani->ScalerVectorMinusVector(x1, betatemp, delta, z, delta);
		e_Pd = betatemp * (e_Pd + alphatemp * d_Pd);
		d_Pd = z_r + betatemp * betatemp * d_Pd;
	}
	innerIter = j;
};

void SolversTR::PreConditioner(Variable *x, Vector *eta, Vector *result)
{
	// default one means no preconditioner.
	eta->CopyTo(result);
};

void SolversTR::PrintGenInfo(void)
{
	Solvers::PrintGenInfo();
	printf("nH:%d,rho:%.2e,radius:%.3e,tCGstatus:%s,innerIter:%d,", nH, rho, Delta, tCGstatusSetnames[tCGstatus].c_str(), innerIter);
};

void SolversTR::CheckParams(void)
{
	Solvers::CheckParams();

	char YES[] = "YES";
	char NO[] = "NO";
	char *status;

	std::cout << "TRUST REGION TYPE METHODS PARAMETERS:" << std::endl;
	status = (initial_Delta > 0) ? YES : NO;
	std::cout << "initial_Delta :" << std::setw(15) << initial_Delta << "[" << status << "],\t";
	status = (Acceptence_Rho > 0 && Acceptence_Rho < 0.25) ? YES : NO;
	std::cout << "Acceptence_Rho:" << std::setw(15) << Acceptence_Rho << "[" << status << "]" << std::endl;
	status = (Shrinked_tau > 0 && Shrinked_tau < 1) ? YES : NO;
	std::cout << "Shrinked_tau  :" << std::setw(15) << Shrinked_tau << "[" << status << "],\t";
	status = (Magnified_tau > 1) ? YES : NO;
	std::cout << "Magnified tau :" << std::setw(15) << Magnified_tau << "[" << status << "]" << std::endl;
	status = (minimum_Delta > 0 && minimum_Delta <= maximum_Delta) ? YES : NO;
	std::cout << "minimum_Delta :" << std::setw(15) << minimum_Delta << "[" << status << "],\t";
	status = (maximum_Delta > 0 && maximum_Delta >= minimum_Delta) ? YES : NO;
	std::cout << "maximum_Delta :" << std::setw(15) << maximum_Delta << "[" << status << "]" << std::endl;
	status = (Min_Inner_Iter >= 0 && Min_Inner_Iter <= Max_Inner_Iter) ? YES : NO;
	std::cout << "Min_Inner_Iter:" << std::setw(15) << Min_Inner_Iter << "[" << status << "],\t";
	status = (Max_Inner_Iter >= 0 && Max_Inner_Iter >= Min_Inner_Iter) ? YES : NO;
	std::cout << "Max_Inner_Iter:" << std::setw(15) << Max_Inner_Iter << "[" << status << "]" << std::endl;
	status = (theta >= 0) ? YES : NO;
	std::cout << "theta         :" << std::setw(15) << theta << "[" << status << "],\t";
	status = (kappa > 0 && kappa < 1) ? YES : NO;
	std::cout << "kappa         :" << std::setw(15) << kappa << "[" << status << "]" << std::endl;
	status = YES;
	std::cout << "useRand       :" << std::setw(15) << useRand << "[" << status << "]" << std::endl;
};

void SolversTR::UpdateData(void)
{
};

void SolversTR::Acceptence(void)
{
	Prob->Grad(x2, gf2); ng++;
};

void SolversTR::Initialization(const Problem *prob, const Variable *initialx, const Vector *EMPTYETA)
{
	nH = 0;
	Solvers::Initialization(prob, initialx, EMPTYETA);
	eta1 = EMPTYETA->ConstructEmpty();
	eta2 = EMPTYETA->ConstructEmpty();
	zeta = EMPTYETA->ConstructEmpty();
	r = EMPTYETA->ConstructEmpty();
	z = EMPTYETA->ConstructEmpty();
	delta = EMPTYETA->ConstructEmpty();
	Hd = EMPTYETA->ConstructEmpty();

	Acceptence_Rho = 0.1;
	Shrinked_tau = 0.25;
	Magnified_tau = 2;
	minimum_Delta = std::numeric_limits<double>::epsilon();
	maximum_Delta = 1000;
	useRand = false;
	Max_Inner_Iter = 1000;
	Min_Inner_Iter = 0;
	theta = 1;
	kappa = 0.1;
	initial_Delta = 1;
	tCGstatusSetnames = new std::string [TCGSTATUSSETLENGTH];
	tCGstatusSetnames[TR_NEGCURVTURE].assign("NEGCURVTURE");
	tCGstatusSetnames[TR_EXCREGION].assign("EXCREGION");
	tCGstatusSetnames[TR_LCON].assign("LCON");
	tCGstatusSetnames[TR_SCON].assign("SCON");
	tCGstatusSetnames[TR_MAXITER].assign("MAXITER");
};

SolversTR::~SolversTR(void)
{
	delete eta1;
	delete eta2;
	delete zeta;
	delete r;
	delete z;
	delete delta;
	delete Hd;

	delete[] tCGstatusSetnames;
};

void SolversTR::SetParams(PARAMSMAP params)
{
	Solvers::SetParams(params);
	PARAMSMAP::iterator iter;
	for (iter = params.begin(); iter != params.end(); iter++)
	{
		if (iter->first == static_cast<std::string> ("Acceptence_Rho"))
		{
			Acceptence_Rho = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("Shrinked_tau"))
		{
			Shrinked_tau = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("Magnified_tau"))
		{
			Magnified_tau = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("minimum_Delta"))
		{
			minimum_Delta = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("maximum_Delta"))
		{
			maximum_Delta = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("useRand"))
		{
			useRand = ((static_cast<integer> (iter->second)) != 0);
		}
		else
		if (iter->first == static_cast<std::string> ("Max_Inner_Iter"))
		{
			Max_Inner_Iter = static_cast<integer> (iter->second);
		}
		else
		if (iter->first == static_cast<std::string> ("Min_Inner_Iter"))
		{
			Min_Inner_Iter = static_cast<integer> (iter->second);
		}
		else
		if (iter->first == static_cast<std::string> ("theta"))
		{
			theta = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("kappa"))
		{
			kappa = iter->second;
		}
		else
		if (iter->first == static_cast<std::string> ("initial_Delta"))
		{
			initial_Delta = iter->second;
		}
	}
};
