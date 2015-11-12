
#include "L2Sphere.h"

L2Sphere::L2Sphere(integer inn)
{
	// public parameter
	HasHHR = false;

	// Only some combinations exist.
	metric = TRAPEZOID;
	retraction = NORMALIZED;
	VecTran = ISOBYHHR;
	IsIntrApproach = false;
	UpdBetaAlone = false;

	// Status of locking condition
	HasLockCon = false;

	// Fixed parameters
	n = inn;
	ExtrinsicDim = n;
	IntrinsicDim = n - 1;
	name.assign("L2Sphere");

	EMPTYEXTR = new L2SphereVector(n);
	EMPTYINTR = new L2SphereVector(IntrinsicDim);
};

L2Sphere::~L2Sphere()
{
	delete EMPTYEXTR;
	delete EMPTYINTR;
};

double L2Sphere::Metric(Variable *x, Vector *etax, Vector *xix) const
{ //Trapezoidal rule
	const double *etaxTV = etax->ObtainReadData();
	const double *xixTV = xix->ObtainReadData();
	integer inc = 1;
	double result = ddot_(&n, const_cast<double *> (etaxTV), &inc, const_cast<double *> (xixTV), &inc);
	result -= etaxTV[0] * xixTV[0] / 2;
	result -= etaxTV[n - 1] * xixTV[n - 1] / 2;
	return result / (n - 1);
};

void L2Sphere::Projection(Variable *x, Vector *v, Vector *result) const
{
	const double *xl = x->ObtainReadData();
	double nume = Metric(x, x, v);
	ScalerVectorAddVector(x, -nume, x, v, result);
};

void L2Sphere::Retraction(Variable *x, Vector *etax, Variable *result) const
{// exponential mapping
	double norm = sqrt(Metric(x, etax, etax));
	if (norm < std::numeric_limits<double>::epsilon())
		ScaleTimesVector(x, cos(norm), x, result);
	else
		VectorLinearCombination(x, cos(norm), x, sin(norm) / norm, etax, result);
};

void L2Sphere::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	xiy->CopyTo(result);
	std::cout << "The cotangent vector has not been implemented!" << std::endl;
};

void L2Sphere::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
{
	if (IsEtaXiSameDir)
	{
		VectorTransport(x, etax, y, xix, result);

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
	std::cout << "Warning: The differentiated retraction has not been implemented!" << std::endl;
	xix->CopyTo(result);
};

double L2Sphere::Beta(Variable *x, Vector *etax) const
{
	return 1;
};

void L2Sphere::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
{
	if (!etax->TempDataExist("xdydn2"))
	{
		Vector *xdy = x->ConstructEmpty();
		SharedSpace *Sharedxdy = new SharedSpace(xdy);
		VectorAddVector(x, x, y, xdy);
		ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
		etax->AddToTempData("xdydn2", Sharedxdy);
	}
	const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
	Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
	ScalerVectorAddVector(x, -2.0 * Metric(x, xix, y), xdydn2, xix, result);
};

void L2Sphere::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	if (!etax->TempDataExist("xdydn2"))
	{
		Vector *xdy = x->ConstructEmpty();
		SharedSpace *Sharedxdy = new SharedSpace(xdy);
		VectorAddVector(x, x, y, xdy);
		ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
		etax->AddToTempData("xdydn2", Sharedxdy);
	}
	const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
	Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
	ScalerVectorAddVector(x, -2.0 * Metric(x, xiy, x), xdydn2, xiy, result);
};

void L2Sphere::HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
{
	if (!etax->TempDataExist("xdydn2"))
	{
		Vector *xdy = x->ConstructEmpty();
		SharedSpace *Sharedxdy = new SharedSpace(xdy);
		VectorAddVector(x, x, y, xdy);
		ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
		etax->AddToTempData("xdydn2", Sharedxdy);
	}
	const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
	Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
	const double *xdydn2TV = xdydn2->ObtainReadData();

	integer ell = Hx->Getsize()[0];
	integer length = etax->Getlength();
	const double *M = Hx->ObtainReadData();
	double *Hxpy = new double[ell];

	char *transn = const_cast<char *> ("n");
	double one = 1, zero = 0;
	integer inc = 1, N = ell;
	dgemv_(transn, &N, &length, &one, const_cast<double *> (M + start * N), &N, const_cast<double *> (xdydn2TV), &inc, &zero, Hxpy, &inc);

	double scaler = -2.0;
	Hx->CopyTo(result);

	Variable *xflat = x->ConstructEmpty();
	x->CopyTo(xflat);
	double *xflatptr = xflat->ObtainWritePartialData();
	xflatptr[0] /= (2 * (n - 1));
	xflatptr[n - 1] /= (2 * (n - 1));
	for (integer i = 1; i < n - 1; i++)
	{
		xflatptr[i] /= (n - 1);
	}
	double *resultL = result->ObtainWritePartialData();
	dger_(&length, &N, &scaler, Hxpy, &inc, xflatptr, &inc, resultL + start * N, &N);
	delete[] Hxpy;
	delete xflat;
};

void L2Sphere::TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
{
	if (!etax->TempDataExist("xdydn2"))
	{
		Vector *xdy = x->ConstructEmpty();
		SharedSpace *Sharedxdy = new SharedSpace(xdy);
		VectorAddVector(x, x, y, xdy);
		ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
		etax->AddToTempData("xdydn2", Sharedxdy);
	}

	integer ell = Hx->Getsize()[0];
	integer length = etax->Getlength();
	const double *M = Hx->ObtainReadData();
	double *Hty = new double[ell];

	Variable *yflat = y->ConstructEmpty();
	y->CopyTo(yflat);
	double *yflatptr = yflat->ObtainWritePartialData();
	yflatptr[0] /= (2 * (n - 1));
	yflatptr[n - 1] /= (2 * (n - 1));
	for (integer i = 1; i < n - 1; i++)
	{
		yflatptr[i] /= (n - 1);
	}

	char *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	integer inc = 1, N = ell;
	dgemv_(transt, &length, &N, &one, const_cast<double *> (M + start), &N, yflatptr, &inc, &zero, Hty, &inc);

	double scaler = -2.0;
	Hx->CopyTo(result);


	const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
	Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
	const double *xdydn2TV = xdydn2->ObtainReadData();

	double *resultL = result->ObtainWritePartialData();
	dger_(&length, &N, &scaler, const_cast<double *> (xdydn2TV), &inc, Hty, &inc, resultL + start, &N);
	delete[] Hty;
	delete yflat;
};

void L2Sphere::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
{
	HInvTran(x, etax, y, Hx, 0, etax->Getlength(), result);
	TranH(x, etax, y, result, 0, etax->Getlength(), result);
};

void L2Sphere::ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const
{
	etax->CopyTo(etaxflat);
	double *etaxflatTV = etaxflat->ObtainWritePartialData();
	double intv = 1.0 / (n - 1);
	ScaleTimesVector(x, intv, etaxflat, etaxflat);
	etaxflatTV[0] /= 2;
	etaxflatTV[n - 1] /= 2;
};

void L2Sphere::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
{
	std::cout << "Routine of obtaining intrinsic representations has not been done!" << std::endl;
};

void L2Sphere::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
{
	std::cout << "Routine of obtaining extrinsic representations has not been done!" << std::endl;
};

void L2Sphere::IntrProjection(Variable *x, Vector *v, Vector *result) const
{
	v->CopyTo(result);
};

void L2Sphere::ExtrProjection(Variable *x, Vector *v, Vector *result) const
{
	const double *xl = x->ObtainReadData();
	double nume = Metric(x, x, v);
	ScalerVectorAddVector(x, -nume, x, v, result);
};

void L2Sphere::CheckParams(void) const
{
	std::string Repa4NSMetricnames[L2SPHEREMETRICLENGTH] = { "TRAPEZOID" };
	std::string Repa4NSRetractionnames[L2SPHERERETRACTIONLENGTH] = { "NORMALIZED" };
	std::string Repa4NSVectorTransportnames[L2SPHEREVECTORTRANSPORTLENGTH] = { "ISOBYHHR" };
	Manifold::CheckParams();
	std::cout << name << " PARAMETERS:" << std::endl;
	std::cout << "n             :" << std::setw(15) << n << ",\t";
	std::cout << "metric        :" << std::setw(15) << Repa4NSMetricnames[metric] << std::endl; 
	std::cout << "retraction    :" << std::setw(15) << Repa4NSRetractionnames[retraction] << ",\t";
	std::cout << "VecTran       :" << std::setw(15) << Repa4NSVectorTransportnames[VecTran] << std::endl;
};

void L2Sphere::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
{
	//egf->CopyTo(gf);//--
	//return;//--

	if (prob->GetUseHess())
	{
		Vector *segf = egf->ConstructEmpty();
		segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
		egf->CopyTo(segf);
		SharedSpace *Sharedegf = new SharedSpace(segf);
		x->AddToTempData("EGrad", Sharedegf);
	}
	ExtrProjection(x, egf, gf);
};

void L2Sphere::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
{
	//exix->CopyTo(xix);//---
	//return; //---
	const double *xptr = x->ObtainReadData();
	Variable *xcubed = x->ConstructEmpty();
	SharedSpace *Sharedxcubed = new SharedSpace(xcubed);
	double *xcubedptr = xcubed->ObtainWriteEntireData();
	for (integer i = 0; i < n; i++)
	{
		xcubedptr[i] = xptr[i] * xptr[i] * xptr[i];
	}
	double a1 = Metric(x, xcubed, xcubed);

	const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
	Vector *egfVec = Sharedegf->GetSharedElement();
	double a2 = Metric(x, egfVec, xcubed);

	Vector *x2etax = etax->ConstructEmpty();
	double *x2etaxptr = x2etax->ObtainWriteEntireData();
	const double *etaxptr = etax->ObtainReadData();
	for (integer i = 0; i < n; i++)
	{
		x2etaxptr[i] = xptr[i] * xptr[i] * etaxptr[i];
	}
	ScalerVectorAddVector(x, -3.0 * a2 / a1, x2etax, exix, xix);
	delete x2etax;
	ExtrProjection(x, xix, xix);
};
