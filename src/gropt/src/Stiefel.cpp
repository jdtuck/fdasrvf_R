
#include "Stiefel.h"

//==============================================================================================
Stiefel::Stiefel(integer inn, integer inp)
{
	// public parameter
	HasHHR = false;

	// Only some combinations exist.
	metric = EUCLIDEAN;
	retraction = QF;
	VecTran = PARALLELIZATION;
	IsIntrApproach = true;
	UpdBetaAlone = false;

	// Status of locking condition
	HasLockCon = false;

	// Fixed parameters
	n = inn;
	p = inp;
	ExtrinsicDim = n * p;
	IntrinsicDim = n * p - p * (p + 1) / 2;
	name.assign("Stiefel");
	EMPTYEXTR = new StieVector(n, p);
	EMPTYINTR = new StieVector(IntrinsicDim);
};

Stiefel::~Stiefel(void)
{
	delete EMPTYEXTR;
	delete EMPTYINTR;
};

// Choose the default parameters
void Stiefel::ChooseStieParamsSet1(void)
{
	metric = EUCLIDEAN;
	retraction = QF;
	VecTran = PARALLELIZATION;
	IsIntrApproach = true;
	UpdBetaAlone = false;
	HasHHR = false;
	HasLockCon = false;
};

// Choose the parameters
void Stiefel::ChooseStieParamsSet2(void)
{
	metric = EUCLIDEAN;
	retraction = CONSTRUCTED;
	VecTran = PARALLELIZATION;
	IsIntrApproach = true;
	UpdBetaAlone = false;
	HasHHR = false;
	HasLockCon = true;
};

void Stiefel::CheckParams(void) const
{
	std::string StieMetricnames[STIEMETRICLENGTH] = { "EUCLIDEAN", "CANONICAL" };
	std::string StieRetractionnames[STIERETRACTIONLENGTH] = { "QF", "POLAR", "EXP", "CONSTRUCTED" };
	std::string StieVectorTransportnames[STIEVECTORTRANSPORTLENGTH] = { "PARALLELIZATION", "RIGGING", "PARALLELTRANSLATION" };
	Manifold::CheckParams();
	std::cout << name << " PARAMETERS:" << std::endl;
	std::cout << "n             :" << std::setw(15) << n << ",\t";
	std::cout << "p             :" << std::setw(15) << p << std::endl;
	std::cout << "metric        :" << std::setw(15) << StieMetricnames[metric] << ",\t";
	std::cout << "retraction    :" << std::setw(15) << StieRetractionnames[retraction] << std::endl;
	std::cout << "VecTran       :" << std::setw(15) << StieVectorTransportnames[VecTran] << std::endl;
};

void Stiefel::IntrProjection(Variable *x, Vector *etax, Vector *result) const
{
	etax->CopyTo(result);
};

void Stiefel::ExtrProjection(Variable *x, Vector *v, Vector *result) const
{
	integer N = n, P = p, inc = 1, Length = N * P;
	double *symUtV = new double[P * P];
	const double *U = x->ObtainReadData();
	const double *V = v->ObtainReadData();
	double *resultTV = result->ObtainWriteEntireData();

	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (U), &N, const_cast<double *> (V), &N, &zero, symUtV, &P);
	for (integer i = 0; i < P; i++)
	{
		for (integer j = i + 1; j < P; j++)
		{
			symUtV[i + j * P] += symUtV[j + i * P];
			symUtV[i + j * P] /= 2.0;
			symUtV[j + i * P] = symUtV[i + j * P];
		}
	}
	if (V != resultTV)
		dcopy_(&Length, const_cast<double *> (V), &inc, resultTV, &inc);
	double negone = -1;
	dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (U), &N, symUtV, &P, &one, resultTV, &N);
	delete [] symUtV;
};

//=================================================================================
double Stiefel::Metric(Variable *x, Vector *etax, Vector *xix) const
{
	if (metric == EUCLIDEAN)
		return Manifold::Metric(x, etax, xix);
	std::cout << "Error: Metric has not been done!" << std::endl;
    exit(0);
	return 0;
};

void Stiefel::Projection(Variable *x, Vector *v, Vector *result) const
{
	if (IsIntrApproach)
		IntrProjection(x, v, result);
	else
		ExtrProjection(x, v, result);
};

void Stiefel::Retraction(Variable *x, Vector *etax, Variable *result) const
{
	if (retraction == QF)
		return qfRetraction(x, etax, result);

	if (retraction == CONSTRUCTED)
		return ConRetraction(x, etax, result);

	std::cout << "Error: Retraction has not been done!" << std::endl;
};

void Stiefel::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	if (retraction == QF)
		return qfcoTangentVector(x, etax, y, xiy, result);
	
	if (retraction == CONSTRUCTED)
		return ConcoTangentVector(x, etax, y, xiy, result);

	std::cout << "Error: coTangentVector has not been done!" << std::endl;
};

void Stiefel::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
{
	if (retraction == QF)
		return DiffqfRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

	if (retraction == CONSTRUCTED)
		return DiffConRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

	std::cout << "Error: DiffRetraction has not been done!" << std::endl;
};

double Stiefel::Beta(Variable *x, Vector *etax) const
{
	if (! HasHHR && ! UpdBetaAlone)
		return 1;

	if (! etax->TempDataExist("beta"))
	{
		Variable *y = x->ConstructEmpty();
		Vector *xiy = etax->ConstructEmpty();
		Retraction(x, etax, y);
		DiffRetraction(x, etax, y, etax, xiy, true);
		delete y;
		delete xiy;
	}

	const SharedSpace *beta = etax->ObtainReadTempData("beta");
	const double *betav = beta->ObtainReadData();
	return betav[0];
};

void Stiefel::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
{
	if (VecTran == PARALLELIZATION && !HasHHR)
		return Manifold::VectorTransport(x, etax, y, xix, result);

	if (HasHHR)
		return LCVectorTransport(x, etax, y, xix, result);

	std::cout << "Error: VectorTransport has not been done!" << std::endl;
};

void Stiefel::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	if (VecTran == PARALLELIZATION && !HasHHR)
		return Manifold::InverseVectorTransport(x, etax, y, xiy, result);

	if (HasHHR)
		return LCInverseVectorTransport(x, etax, y, xiy, result);

	std::cout << "Error: InverseVectorTransport has not been done!" << std::endl;
};

void Stiefel::HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
{
	if (VecTran == PARALLELIZATION && !HasHHR)
		return Manifold::HInvTran(x, etax, y, Hx, start, end, result);

	if (HasHHR)
		return LCHInvTran(x, etax, y, Hx, start, end, result);

	std::cout << "Error: HInvTran has not been done!" << std::endl;
};

void Stiefel::TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
{
	if (VecTran == PARALLELIZATION && !HasHHR)
		return Manifold::TranH(x, etax, y, Hx, start, end, result);

	if (HasHHR)
		return LCTranH(x, etax, y, Hx, start, end, result);

	std::cout << "Error: TranH has not been done!" << std::endl;
};

void Stiefel::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
{
	if (VecTran == PARALLELIZATION && !HasHHR)
		return Manifold::TranHInvTran(x, etax, y, Hx, result);

	if (HasHHR)
		return LCTranHInvTran(x, etax, y, Hx, result);

	std::cout << "Error: TranHInvTran has not been done!" << std::endl;
};

void Stiefel::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
{
	if (metric == EUCLIDEAN)
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		ExtrProjection(x, egf, gf);
		return;
	}
	std::cout << "Warning:The function converting Eucidean Gradient to Riemannian Gradient has not been done!" << std::endl;
};

void Stiefel::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector *xix, const Problem *prob) const
{
	if (metric == EUCLIDEAN)
	{
		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		integer inc = 1, N = n, P = p, Length = N * P;
		double *symxtegfptr;
		SharedSpace *symxtegf;
		if (x->TempDataExist("symxtegf"))
		{
			symxtegf = const_cast<SharedSpace *> (x->ObtainReadTempData("symxtegf"));
			symxtegfptr = const_cast<double *> (symxtegf->ObtainReadData());
		}
		else
		{
			const double *xxM = x->ObtainReadData();
			const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
			Vector *egfVec = Sharedegf->GetSharedElement();
			const double *egf = egfVec->ObtainReadData();
			symxtegf = new SharedSpace(2, p, p);
			symxtegfptr = symxtegf->ObtainWriteEntireData();
			dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xxM), &N, const_cast<double *> (egf), &N, &zero, symxtegfptr, &P);
			for (integer i = 0; i < p; i++)
			{
				for (integer j = i + 1; j < p; j++)
				{
					symxtegfptr[i + j * p] += symxtegfptr[j + i * p];
					symxtegfptr[i + j * p] /= 2.0;
					symxtegfptr[j + i * p] = symxtegfptr[i + j * p];
				}
			}
		}

		exix->CopyTo(xix);
		double *resultTV = xix->ObtainWritePartialData();
		const double *etaxTV = etax->ObtainReadData();

		double negone = -1;
		dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (etaxTV), &N, symxtegfptr, &P, &one, resultTV, &N);
		ExtrProjection(x, xix, xix);
		if (!x->TempDataExist("symxtegf"))
		{
			x->AddToTempData("symxtegf", symxtegf);
		}
		return;
	}
	std::cout << "Warning:The function converting action of Eucidean Hessian to action of Riemannian Hessian has not been done!" << std::endl;
};

void Stiefel::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
{
	if (retraction == QF)
		ObtainIntrHHR(x, etax, result);
	else
	if (retraction == CONSTRUCTED)
		ObtainIntrSquare(x, etax, result);
	else
		std::cout << "Warning: computing intrinsinc representation from extrinsic has not been implemented!" << std::endl;
};

void Stiefel::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
{
	if (retraction == QF)
		ObtainExtrHHR(x, intretax, result);
	else
	if (retraction == CONSTRUCTED)
		ObtainExtrSquare(x, intretax, result);
	else
		std::cout << "Warning: computing extrinsic representation from intrinsinc has not been implemented!" << std::endl;
};

void Stiefel::qfRetraction(Variable *x, Vector *etax, Variable *result) const
{
	//x->Print("x in qf:");//---
	const double *U = x->ObtainReadData();
	const double *V;
	Vector *exetax = nullptr;
	if (IsIntrApproach)
	{
		exetax = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		V = exetax->ObtainReadData();
		//exetax->Print("exetax:");//---
	}
	else
	{
		V = etax->ObtainReadData();
	}
	double *resultM = result->ObtainWriteEntireData();
	SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
	double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
	SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
	double *tau = HHRTau->ObtainWriteEntireData();

	integer N = x->Getsize()[0], P = x->Getsize()[1], Length = N * P, inc = 1;
	double one = 1, zero = 0;
	dcopy_(&Length, const_cast<double *> (V), &inc, ptrHHR, &inc);
	daxpy_(&Length, &one, const_cast<double *> (U), &inc, ptrHHR, &inc);
	integer *jpvt = new integer[P];
	integer info;
	integer lwork = -1;
	double lworkopt;
	for (integer i = 0; i < P; i++)
		jpvt[i] = i + 1;
	dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
	lwork = static_cast<integer> (lworkopt);
	double *work = new double[lwork];
	dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
	if (info < 0)
		std::cout << "Error in qr decomposition!" << std::endl;
	for (integer i = 0; i < P; i++)
	{
		if (jpvt[i] != (i + 1))
			std::cout << "Error in qf retraction!" << std::endl;
	}
	double *signs = new double[P];
	for (integer i = 0; i < P; i++)
		signs[i] = (ptrHHR[i + i * N] >= 0) ? 1 : -1;
	dcopy_(&Length, ptrHHR, &inc, resultM, &inc);
	dorgqr_(&N, &P, &P, resultM, &N, tau, work, &lwork, &info);
	if (info < 0)
		std::cout << "Error in forming Q matrix!" << std::endl;
	for (integer i = 0; i < P; i++)
		dscal_(&N, signs + i, resultM + i * N, &inc);
	result->AddToTempData("HHR", HouseHolderResult);
	result->AddToTempData("HHRTau", HHRTau);
	delete[] jpvt;
	delete[] work;
	delete[] signs;
	if (exetax != nullptr)
		delete exetax;
	//result->Print("result in qf:");//---
};

void Stiefel::qfcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	const double *yM = y->ObtainReadData();
	Vector *exresult = EMPTYEXTR->ConstructEmpty();
	double *exresultTV = exresult->ObtainWriteEntireData();

	Vector *extempy = nullptr;
	const double *extempyTV;
	if (IsIntrApproach)
	{
		extempy = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(y, xiy, extempy);
		extempyTV = extempy->ObtainReadData();
	}
	else
	{
		extempyTV = xiy->ObtainReadData();
	}
	double *ytxiy = new double[p * p];

	char *transt = const_cast<char *> ("t"), *transn = const_cast<char *> ("n");
	integer N = n, P = p, inc = 1;
	double one = 1, zero = 0;
	dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (yM), &N, const_cast<double *> (extempyTV), &N, &zero, ytxiy, &P);
	for (integer i = 0; i < p; i++)
	{
		for (integer j = i; j < p; j++)
		{
			ytxiy[i + j * p] = -ytxiy[i + j * p];
		}
	}
	dgemm_(transn, transn, &N, &P, &P, &one, const_cast<double *> (yM), &N, ytxiy, &P, &zero, exresultTV, &N);
	integer Length = N * P;
	daxpy_(&Length, &one, const_cast<double *> (extempyTV), &inc, exresultTV, &inc);

	const SharedSpace *HHR = y->ObtainReadTempData("HHR");
	const double *ptrHHR = HHR->ObtainReadData();
	double sign;
	for (integer i = 0; i < P; i++)
	{
		sign = (ptrHHR[i + i * N] >= 0) ? 1 : -1;
		dscal_(&N, &sign, exresultTV + i * N, &inc);
	}

	char *left = const_cast<char *> ("r"), *up = const_cast<char *> ("u"), *nonunit = const_cast<char *> ("n");
	dtrsm_(left, up, transt, nonunit, &N, &P, &one, const_cast<double *> (ptrHHR), &N, exresultTV, &N);

	ExtrProjection(x, exresult, exresult);
	if (IsIntrApproach)
	{
		ObtainIntr(x, exresult, result);
	}
	else
	{
		exresult->CopyTo(result);
	}

	delete[] ytxiy;
	delete exresult;
	if (extempy != nullptr)
		delete extempy;
};

void Stiefel::DiffqfRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
{
	Vector *extempx = EMPTYEXTR->ConstructEmpty();
	const double *extempxTV;
	if (IsIntrApproach)
	{
		ObtainExtr(x, xix, extempx);
		extempxTV = extempx->ObtainReadData();
	}
	else
	{
		xix->CopyTo(extempx);
		extempxTV = extempx->ObtainWritePartialData();
	}
	const double *yM = y->ObtainReadData();
	double *resultTV = result->ObtainWriteEntireData();
	const SharedSpace *HHR = y->ObtainReadTempData("HHR");
	const double *ptrHHR = HHR->ObtainReadData();
	double *YtVRinv = new double[p * p];
	integer inc = 1, N = n, P = p;
	char *left = const_cast<char *> ("r"), *up = const_cast<char *> ("u"),
		*transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t"), *nonunit = const_cast<char *> ("n");
	double one = 1, zero = 0;
	dtrsm_(left, up, transn, nonunit, &N, &P, &one, const_cast<double *> (ptrHHR), &N, const_cast<double *> (extempxTV), &N);

	double sign;
	for (integer i = 0; i < P; i++)
	{
		sign = (ptrHHR[i + i * N] >= 0) ? 1 : -1;
		dscal_(&N, &sign, const_cast<double *> (extempxTV + i * N), &inc);
	}
	dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (yM), &N, const_cast<double *> (extempxTV), &N, &zero, YtVRinv, &P);
	for (integer i = 0; i < p; i++)
	{
		YtVRinv[i + p * i] = -YtVRinv[i + p * i];
		for (integer j = i + 1; j < p; j++)
		{
			YtVRinv[i + p * j] = -YtVRinv[j + p * i] - YtVRinv[i + p * j];
			YtVRinv[j + p * i] = 0;
		}
	}
	dgemm_(transn, transn, &N, &P, &P, &one, const_cast<double *> (yM), &N, YtVRinv, &P, &one, const_cast<double *> (extempxTV), &N);
	if (IsIntrApproach)
	{
		ObtainIntr(y, extempx, result);
	}
	else
	{
		extempx->CopyTo(result);
	}
	delete[] YtVRinv;
	delete extempx;

	if (IsEtaXiSameDir && (HasHHR || UpdBetaAlone))
	{
		const double *etaxTV = etax->ObtainReadData();
		const double *xixTV = xix->ObtainReadData();
		double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
		SharedSpace *beta = new SharedSpace(1, 3);
		double *betav = beta->ObtainWriteEntireData();
		betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
		betav[1] = Metric(x, etax, etax);
		betav[2] = Metric(x, result, result) * EtatoXi * EtatoXi;
		etax->AddToTempData("beta", beta);

		if (HasHHR)
		{
			Vector *TReta = result->ConstructEmpty();
			result->CopyTo(TReta);
			ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
			SharedSpace *SharedTReta = new SharedSpace(TReta);
			etax->AddToTempData("betaTReta", SharedTReta);
		}
	}
};

void Stiefel::ObtainIntrHHR(Variable *x, Vector *etax, Vector *result) const
{
	if (!x->TempDataExist("HHR"))
	{
		const double *xM = x->ObtainReadData();
		SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
		double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
		SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
		double *tau = HHRTau->ObtainWriteEntireData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], Length = N * P, inc = 1;
		double one = 1, zero = 0;
		dcopy_(&Length, const_cast<double *> (xM), &inc, ptrHHR, &inc);
		integer *jpvt = new integer[P];
		integer info;
		integer lwork = -1;
		double lworkopt;
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		for (integer i = 0; i < P; i++)
			jpvt[i] = i + 1;
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
		x->AddToTempData("HHR", HouseHolderResult);
		x->AddToTempData("HHRTau", HHRTau);
		if (info < 0)
			std::cout << "Error in qr decomposition!" << std::endl;
		for (integer i = 0; i < P; i++)
		{
			if (jpvt[i] != (i + 1))
				std::cout << "Error in qf retraction!" << std::endl;
		}
		delete[] jpvt;
		delete[] work;
	}

	const double *xM = x->ObtainReadData();
	const double *etaxTV = etax->ObtainReadData();
	const SharedSpace *HHR = x->ObtainReadTempData("HHR");
	const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
	double *resultTV = result->ObtainWriteEntireData();
	const double *ptrHHR = HHR->ObtainReadData();
	const double *ptrHHRTau = HHRTau->ObtainReadData();

	char *transt = const_cast<char *> ("t"), *sidel = const_cast<char *> ("l");
	integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
	integer info;
	integer lwork = -1;
	double lworkopt;
	double *tempspace = new double[n * p];
	dormqr_(sidel, transt, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), tempspace, &N, &lworkopt, &lwork, &info);
	lwork = static_cast<integer> (lworkopt);
	double *work = new double[lwork];
	dcopy_(&Length, const_cast<double *> (etaxTV), &inc, tempspace, &inc);
	dormqr_(sidel, transt, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), tempspace, &N, work, &lwork, &info);
	double sign;
	for (integer i = 0; i < p; i++)
	{
		sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
		dscal_(&P, &sign, tempspace + i, &N);
	}
	double r2 = sqrt(2.0);
	integer idx = 0;
	for (integer i = 0; i < p; i++)
	{
		for (integer j = i + 1; j < p; j++)
		{
			resultTV[idx] = r2 * tempspace[j + i * n];
			idx++;
		}
	}

	for (integer i = 0; i < p; i++)
	{
		for (integer j = p; j < n; j++)
		{
			resultTV[idx] = tempspace[j + i * n];
			idx++;
		}
	}

	delete[] work;
	delete[] tempspace;
};

void Stiefel::ObtainExtrHHR(Variable *x, Vector *intretax, Vector *result) const
{
	if (!x->TempDataExist("HHR"))
	{
		const double *xM = x->ObtainReadData();
		SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
		double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
		SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
		double *tau = HHRTau->ObtainWriteEntireData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], Length = N * P, inc = 1;
		double one = 1, zero = 0;
		dcopy_(&Length, const_cast<double *> (xM), &inc, ptrHHR, &inc);
		integer *jpvt = new integer[P];
		integer info;
		integer lwork = -1;
		double lworkopt;
		for (integer i = 0; i < P; i++)
			jpvt[i] = i + 1;
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
		x->AddToTempData("HHR", HouseHolderResult);
		x->AddToTempData("HHRTau", HHRTau);
		if (info < 0)
			std::cout << "Error in qr decomposition!" << std::endl;
		for (integer i = 0; i < P; i++)
		{
			if (jpvt[i] != (i + 1))
				std::cout << "Error in qf retraction!" << std::endl;
		}
		delete[] jpvt;
		delete[] work;
	}

	const double *xM = x->ObtainReadData();
	const SharedSpace *HHR = x->ObtainReadTempData("HHR");
	const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
	const double *ptrHHR = HHR->ObtainReadData();
	const double *ptrHHRTau = HHRTau->ObtainReadData();
	const double *intretaxTV = intretax->ObtainReadData();
	double *resultTV = result->ObtainWriteEntireData();

	char *transn = const_cast<char *> ("n"), *sidel = const_cast<char *> ("l");
	integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
	integer info;
	double r2 = sqrt(2.0);
	integer idx = 0;
	for (integer i = 0; i < p; i++)
	{
		resultTV[i + i * n] = 0;
		for (integer j = i + 1; j < p; j++)
		{
			resultTV[j + i * n] = intretaxTV[idx] / r2;
			resultTV[i + j * n] = -resultTV[j + i * n];
			idx++;
		}
	}

	for (integer i = 0; i < p; i++)
	{
		for (integer j = p; j < n; j++)
		{
			resultTV[j + i * n] = intretaxTV[idx];
			idx++;
		}
	}

	double sign;
	for (integer i = 0; i < p; i++)
	{
		sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
		dscal_(&P, &sign, resultTV + i, &N);
	}
	integer lwork = -1;
	double lworkopt;
	dormqr_(sidel, transn, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), resultTV, &N, &lworkopt, &lwork, &info);
	lwork = static_cast<integer> (lworkopt);

	double *work = new double[lwork];
	dormqr_(sidel, transn, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), resultTV, &N, work, &lwork, &info);
	delete[] work;
};

void Stiefel::ConRetraction(Variable *x, Vector *etax, Variable *result) const
{ // only accept intrinsic approach
	const double *V = etax->ObtainReadData();
	Vector *exetax = nullptr;
	double *M = new double[3 * n * n + 2 * n];
	double *wr = M + n * n;
	double *wi = wr + n;
	double *Vs = wi + n;
	double *VsT = Vs + n * n;

	double r2 = sqrt(2.0);
	integer idx = 0;
	for (integer i = 0; i < p; i++)
	{
		M[i + i * n] = 0;
		for (integer j = i + 1; j < p; j++)
		{
			M[j + i * n] = V[idx] / r2;
			M[i + j * n] = -M[j + i * n];
			idx++;
		}
	}

	for (integer i = 0; i < p; i++)
	{
		for (integer j = p; j < n; j++)
		{
			M[j + i * n] = V[idx];
			M[i + j * n] = -V[idx];
			idx++;
		}
	}
	for (integer i = p; i < n; i++)
	{
		for (integer j = p; j < n; j++)
		{
			M[j + i * n] = 0;
		}
	}
	char *jobv = const_cast<char *> ("V"), *sortn = const_cast<char *> ("N");
	integer N = n, P = p, NmP = n - p, sdim, info;
	integer lwork = -1;
	double lworkopt;

	dgees_(jobv, sortn, nullptr, &N, M, &N, &sdim, wr, wi, Vs, &N, &lworkopt, &lwork, nullptr, &info);
	lwork = static_cast<integer> (lworkopt);
	double *work = new double[lwork];
	dgees_(jobv, sortn, nullptr, &N, M, &N, &sdim, wr, wi, Vs, &N, work, &lwork, nullptr, &info);

	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double cosv, sinv;
	integer two = 2, inc = 1;
	double one = 1, zero = 0;
	double block[4];
	for (integer i = 0; i < n; i++)
	{
		if (i + 1 < n && fabs(M[i + (i + 1) * n]) > std::numeric_limits<double>::epsilon())
		{
			cosv = cos(M[i + (i + 1) * n]);
			sinv = sin(M[i + (i + 1) * n]);
			block[0] = cosv; block[1] = -sinv; block[2] = sinv; block[3] = cosv;
			dgemm_(transn, transn, &N, &two, &two, &one, Vs + i * n, &N, block, &two, &zero, VsT + i * n, &N);
			i++;
		}
		else
		{
			dcopy_(&N, Vs + i * n, &inc, VsT + i * n, &inc);
		}
	}
	dgemm_(transn, transt, &N, &N, &N, &one, VsT, &N, Vs, &N, &zero, M, &N);

	if (!x->TempDataExist("Perp"))
	{
		ObtainPerp(x);
	}
	const SharedSpace *SharedSpacePerp = x->ObtainReadTempData("Perp");
	const double *Perp = SharedSpacePerp->ObtainReadData();

	const double *xM = x->ObtainReadData();
	double *resultM = result->ObtainWriteEntireData();
	SharedSpace *ResultSharedPerp = new SharedSpace(2, n, n - p);
	double *ResultPerp = ResultSharedPerp->ObtainWriteEntireData();

	dgemm_(transn, transn, &P, &P, &P, &one, const_cast<double *> (xM), &N, M, &N, &zero, resultM, &N);
	dgemm_(transn, transn, &P, &P, &NmP, &one, const_cast<double *> (Perp), &N, M + p, &N, &one, resultM, &N);

	dgemm_(transn, transn, &NmP, &P, &P, &one, const_cast<double *> (xM + p), &N, M, &N, &zero, resultM + p, &N);
	dgemm_(transn, transn, &NmP, &P, &NmP, &one, const_cast<double *> (Perp + p), &N, M + p, &N, &one, resultM + p, &N);

	dgemm_(transn, transn, &P, &NmP, &P, &one, const_cast<double *> (xM), &N, M + n * p, &N, &zero, ResultPerp, &N);
	dgemm_(transn, transn, &P, &NmP, &NmP, &one, const_cast<double *> (Perp), &N, M + n * p + p, &N, &one, ResultPerp, &N);

	dgemm_(transn, transn, &NmP, &NmP, &P, &one, const_cast<double *> (xM + p), &N, M + n * p, &N, &zero, ResultPerp + p, &N);
	dgemm_(transn, transn, &NmP, &NmP, &NmP, &one, const_cast<double *> (Perp + p), &N, M + n * p + p, &N, &one, ResultPerp + p, &N);

	result->AddToTempData("Perp", ResultSharedPerp);

	delete[] work;
	delete[] M;
};

void Stiefel::ConcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	xiy->CopyTo(result);
	std::cout << "The cotangent vector for the constructed retraction has not been implemented!" << std::endl;
};

void Stiefel::DiffConRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
{
	if (IsEtaXiSameDir)
	{
		Manifold::VectorTransport(x, etax, y, xix, result);

		if (IsEtaXiSameDir && (HasHHR || UpdBetaAlone))
		{
			const double *etaxTV = etax->ObtainReadData();
			const double *xixTV = xix->ObtainReadData();
			double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
			SharedSpace *beta = new SharedSpace(1, 3);
			double *betav = beta->ObtainWriteEntireData();
			betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
			betav[1] = Metric(x, etax, etax);
			betav[2] = Metric(x, result, result) * EtatoXi * EtatoXi;
			etax->AddToTempData("beta", beta);

			if (HasHHR)
			{
				Vector *TReta = result->ConstructEmpty();
				result->CopyTo(TReta);
				ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
				SharedSpace *SharedTReta = new SharedSpace(TReta);
				etax->AddToTempData("betaTReta", SharedTReta);
			}
		}
		return;
	}
	std::cout << "Warning: The differentiated retraction of the constructed retraction has not been implemented!" << std::endl;
	xix->CopyTo(result);
};

void Stiefel::ObtainPerp(Variable *x) const
{
	const double *xM = x->ObtainReadData();
	SharedSpace *SharedSpacePerp = new SharedSpace(2, n, n - p);
	double *Perp = SharedSpacePerp->ObtainWriteEntireData();
	for (integer i = 0; i < n * (n - p); i++)
		Perp[i] = genrand_gaussian();
	double *temp = new double[p * (n - p)];
	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double one = 1, zero = 0, neg_one = -1;
	integer P = p, N = n, NmP = n - p;
	// temp = X^T * Perp
	dgemm_(transt, transn, &P, &NmP, &N, &one, const_cast<double *> (xM), &N, Perp, &N, &zero, temp, &P);
	// Perp = Perp - X * temp
	dgemm_(transn, transn, &N, &NmP, &P, &neg_one, const_cast<double *> (xM), &N, temp, &P, &one, Perp, &N);
	delete[] temp;

	integer *jpvt = new integer[NmP];
	integer lwork = 2 * NmP + (1 + NmP) * INITIALBLOCKSIZE, info;
	double *tau = new double[NmP + lwork];
	double *work = tau + NmP;
	for (integer i = 0; i < NmP; i++)
		jpvt[i] = 0;
	dgeqp3_(&N, &NmP, Perp, &N, jpvt, tau, work, &lwork, &info);
	if (info < 0)
		std::cout << "Error in qr decomposition!" << std::endl;
	dorgqr_(&N, &NmP, &NmP, Perp, &N, tau, work, &lwork, &info);
	if (info < 0)
		std::cout << "Error in forming Q matrix!" << std::endl;
	delete[] jpvt;
	delete[] tau;

	x->AddToTempData("Perp", SharedSpacePerp);
};

void Stiefel::ObtainIntrSquare(Variable *x, Vector *etax, Vector *result) const
{
	if (! x->TempDataExist("Perp"))
	{
		ObtainPerp(x);
	}
	const SharedSpace *SharedSpacePerp = x->ObtainReadTempData("Perp");
	const double *Perp = SharedSpacePerp->ObtainReadData();

	const double *xM = x->ObtainReadData();
	const double *etaxTV = etax->ObtainReadData();
	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	integer N = n, P = p, NmP = n - p;
	double one = 1, zero = 0;
	double *tempspace = new double[n * p];
	dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xM), &N, const_cast<double *> (etaxTV), &N, &zero, tempspace, &N);
	dgemm_(transt, transn, &NmP, &P, &N, &one, const_cast<double *> (Perp), &N, const_cast<double *> (etaxTV), &N, &zero, tempspace + p, &N);

	double *resultTV = result->ObtainWriteEntireData();
	double r2 = sqrt(2.0);
	integer idx = 0;
	for (integer i = 0; i < p; i++)
	{
		for (integer j = i + 1; j < p; j++)
		{
			resultTV[idx] = r2 * tempspace[j + i * n];
			idx++;
		}
	}

	for (integer i = 0; i < p; i++)
	{
		for (integer j = p; j < n; j++)
		{
			resultTV[idx] = tempspace[j + i * n];
			idx++;
		}
	}

	delete[] tempspace;
};

void Stiefel::ObtainExtrSquare(Variable *x, Vector *intretax, Vector *result) const
{
	if (!x->TempDataExist("Perp"))
	{
		ObtainPerp(x);
	}
	const SharedSpace *SharedSpacePerp = x->ObtainReadTempData("Perp");
	const double *Perp = SharedSpacePerp->ObtainReadData();
	const double *intretaxTV = intretax->ObtainReadData();
	double *tempspace = new double[n * p];
	double r2 = sqrt(2.0);
	integer idx = 0;
	for (integer i = 0; i < p; i++)
	{
		tempspace[i + i * n] = 0;
		for (integer j = i + 1; j < p; j++)
		{
			tempspace[j + i * n] = intretaxTV[idx] / r2;
			tempspace[i + j * n] = -tempspace[j + i * n];
			idx++;
		}
	}

	for (integer i = 0; i < p; i++)
	{
		for (integer j = p; j < n; j++)
		{
			tempspace[j + i * n] = intretaxTV[idx];
			idx++;
		}
	}
	double *resultTV = result->ObtainWriteEntireData();
	const double *xM = x->ObtainReadData();

	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	integer N = n, P = p, NmP = n - p;
	double one = 1, zero = 0;

	dgemm_(transn, transn, &N, &P, &P, &one, const_cast<double *> (xM), &N, tempspace, &N, &zero, resultTV, &N);
	dgemm_(transn, transn, &N, &P, &NmP, &one, const_cast<double *> (Perp), &N, tempspace + p, &N, &one, resultTV, &N);

	delete[] tempspace;
};

void Stiefel::SetParams(PARAMSMAP params)
{
	Manifold::SetParams(params);
	PARAMSMAP::iterator iter;
	for (iter = params.begin(); iter != params.end(); iter++)
	{
		if (iter->first == static_cast<std::string> ("ParamSet"))
		{
			switch (static_cast<integer> (iter->second))
			{
				case 1:
					ChooseStieParamsSet1();
					break;
				case 2:
					ChooseStieParamsSet2();
					break;
				default:
					break;
			}
		}
	}
};
