
#include "LowRank.h"

LowRank::LowRank(integer inm, integer inn, integer inr) : ProductManifold(3, 
	new Stiefel(inm, inr), static_cast<integer> (1), new Euclidean(inr, inr), static_cast<integer> (1), new Stiefel(inn, inr), static_cast<integer> (1))
{
	m = inm;
	n = inn;
	r = inr;
	name.assign("LowRank");
	delete EMPTYEXTR;
	delete EMPTYINTR;
	EMPTYEXTR = new LowRankVector(m, r, r, n, r);
	EMPTYINTR = new LowRankVector(m * r - r * (r + 1) / 2, 1, r, n * r - r * (r + 1) / 2, 1);
};

LowRank::~LowRank()
{
	for (integer i = 0; i < numofmani; i++)
	{
		delete manifolds[i];
	}
};

void LowRank::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
{
	LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
	LowRankVector *LRetax = dynamic_cast<LowRankVector *> (etax);
	LowRankVector *LRresult = dynamic_cast<LowRankVector *> (result);
	LRresult->NewMemoryOnWrite();

	manifolds[0]->ObtainIntr(LRx->GetElement(0), LRetax->GetElement(0), LRresult->GetElement(0));
	manifolds[1]->ObtainIntr(LRx->GetElement(1), LRetax->GetElement(1), LRresult->GetElement(1));
	manifolds[2]->ObtainIntr(LRx->GetElement(2), LRetax->GetElement(2), LRresult->GetElement(2));

	const double *D = LRx->GetElement(1)->ObtainReadData();
	double *UK = LRresult->GetElement(0)->ObtainWritePartialData() + r * (r - 1) / 2;
	double *VK = LRresult->GetElement(2)->ObtainWritePartialData() + r * (r - 1) / 2;
	double *UOmegaH = LRresult->GetElement(0)->ObtainWritePartialData();
	double *VOmegaH = LRresult->GetElement(2)->ObtainWritePartialData();

	double *UKD = new double[(m - r) * r + (n - r) * r + 2 * r * r];
	double *VKD = UKD + (m - r) * r;
	double *UOmega = VKD + (n - r) * r;
	double *VOmega = UOmega + r * r;
	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	integer MmR = m - r, NmR = n - r, R = r, RR = r * r, inc = 1, MmRR = (m - r) * r, NmRR = (n - r) * r;
	double one = 1, zero = 0;
	dgemm_(transn, transn, &MmR, &R, &R, &one, UK, &MmR, const_cast<double *> (D), &R, &zero, UKD, &MmR);
	dgemm_(transn, transt, &NmR, &R, &R, &one, VK, &NmR, const_cast<double *> (D), &R, &zero, VKD, &NmR);
	dcopy_(&MmRR, UKD, &inc, UK, &inc);
	dcopy_(&NmRR, VKD, &inc, VK, &inc);

	double r2 = sqrt(static_cast<double> (2));
	integer idx = 0;
	for (integer i = 0; i < r; i++)
	{
		UOmega[i + i * r] = 0;
		for (integer j = i + 1; j < r; j++)
		{
			UOmega[j + i * r] = UOmegaH[idx] / r2;
			UOmega[i + j * r] = -UOmega[j + i * r];
			idx++;
		}
	}
	idx = 0;
	for (integer i = 0; i < r; i++)
	{
		VOmega[i + i * r] = 0;
		for (integer j = i + 1; j < r; j++)
		{
			VOmega[j + i * r] = VOmegaH[idx] / r2;
			VOmega[i + j * r] = -VOmega[j + i * r];
			idx++;
		}
	}
	for (integer i = 0; i < r * (r - 1) / 2; i++)
	{
		UOmegaH[i] = 0;
		VOmegaH[i] = 0;
	}
	double *dD = LRresult->GetElement(1)->ObtainWritePartialData();
	dgemm_(transn, transn, &R, &R, &R, &one, UOmega, &R, const_cast<double *> (D), &R, &one, dD, &R);
	dgemm_(transn, transt, &R, &R, &R, &one, const_cast<double *> (D), &R, VOmega, &R, &one, dD, &R);

	delete[] UKD;
};

void LowRank::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
{
	LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
	LowRankVector *LRintretax = dynamic_cast<LowRankVector *> (intretax);
	LowRankVector *LRresult = dynamic_cast<LowRankVector *> (result);
	LRresult->NewMemoryOnWrite();
	LowRankVector *LRintretaxDinv = LRintretax->ConstructEmpty();
	LRintretaxDinv->NewMemoryOnWrite();
	LRintretax->CopyTo(LRintretaxDinv);
	const double *D = LRx->GetElement(1)->ObtainReadData();
	double *UKD = LRintretaxDinv->GetElement(0)->ObtainWritePartialData() + r * (r - 1) / 2;
	double *VKD = LRintretaxDinv->GetElement(2)->ObtainWritePartialData() + r * (r - 1) / 2;
	double *UK = new double[(m - r) * r + (n - r) * r + r * r];
	double *VK = UK + (m - r) * r;
	double *Dinv = VK + (n - r) * r;
	integer *IPIV = new integer[r];

	//- It is not a good approach -------
	integer MmR = m - r, NmR = n - r, R = r, inc = 1, RR = r * r, MmRR = (m - r) * r, NmRR = (n - r) * r, info;
	dcopy_(&RR, const_cast<double *> (D), &inc, Dinv, &inc);
	dgetrf_(&R, &R, Dinv, &R, IPIV, &info);
	integer lwork = -1;
	double lworkopt;
	dgetri_(&R, Dinv, &R, IPIV, &lworkopt, &lwork, &info);
	lwork = static_cast<integer> (lworkopt);
	double *work = new double[lwork];
	dgetri_(&R, Dinv, &R, IPIV, work, &lwork, &info);
	delete[] work;
	delete[] IPIV;
	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	dgemm_(transn, transn, &MmR, &R, &R, &one, UKD, &MmR, Dinv, &R, &zero, UK, &MmR);
	dgemm_(transn, transt, &NmR, &R, &R, &one, VKD, &NmR, Dinv, &R, &zero, VK, &NmR);
	dcopy_(&MmRR, UK, &inc, UKD, &inc);
	dcopy_(&NmRR, VK, &inc, VKD, &inc);

	manifolds[0]->ObtainExtr(LRx->GetElement(0), LRintretaxDinv->GetElement(0), LRresult->GetElement(0));
	manifolds[1]->ObtainExtr(LRx->GetElement(1), LRintretaxDinv->GetElement(1), LRresult->GetElement(1));
	manifolds[2]->ObtainExtr(LRx->GetElement(2), LRintretaxDinv->GetElement(2), LRresult->GetElement(2));

	delete[] UK;
	delete LRintretaxDinv;
};

void LowRank::Retraction(Variable *x, Vector *etax, Variable *result) const
{
	Vector *exetax = EMPTYEXTR->ConstructEmpty();
	ObtainExtr(x, etax, exetax);

	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i]->SetIsIntrApproach(false);
	}
	ProductManifold::Retraction(x, exetax, result);
	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i]->SetIsIntrApproach(true);
	}
	delete exetax;
};

void LowRank::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	Vector *exetax = EMPTYEXTR->ConstructEmpty();
	Vector *exxiy = EMPTYEXTR->ConstructEmpty();
	Vector *exresult = EMPTYEXTR->ConstructEmpty();
	ObtainExtr(x, etax, exetax);
	ObtainExtr(y, xiy, exxiy);

	LowRankVariable *LRy = dynamic_cast<LowRankVariable *> (y);
	LowRankVector *LRexxiy = dynamic_cast<LowRankVector *> (exxiy);

	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i]->SetIsIntrApproach(false);
	}
	double *exxiyU = LRexxiy->GetElement(0)->ObtainWritePartialData();
	double *exxiyD = LRexxiy->GetElement(1)->ObtainWritePartialData();
	double *exxiyV = LRexxiy->GetElement(2)->ObtainWritePartialData();
	const double *Uy = LRy->GetElement(0)->ObtainReadData();
	const double *Dy = LRy->GetElement(1)->ObtainReadData();
	const double *Vy = LRy->GetElement(2)->ObtainReadData();

	double *UxiDDyt = new double[2 * m * r];
	double *Utemp = UxiDDyt + m * r;

	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	integer M = m, R = r, N = n, MR = M * R, NR = N * R, inc = 1;
	double one = 1, zero = 0, negone = -1;
	dgemm_(transn, transn, &M, &R, &R, &one, exxiyU, &M, const_cast<double *> (Dy), &R, &zero, Utemp, &M);
	dgemm_(transn, transt, &M, &R, &R, &one, Utemp, &M, const_cast<double *> (Dy), &R, &zero, exxiyU, &M);
	dgemm_(transn, transn, &M, &R, &R, &one, const_cast<double *> (Uy), &M, exxiyD, &R, &zero, Utemp, &M);
	dgemm_(transn, transt, &M, &R, &R, &one, Utemp, &M, const_cast<double *> (Dy), &R, &zero, UxiDDyt, &M);
	daxpy_(&MR, &one, UxiDDyt, &inc, exxiyU, &inc);

	delete[] UxiDDyt;

	double *VxiDtDy = new double[2 * n * r];
	double *Vtemp = VxiDtDy + n * r;
	dgemm_(transn, transt, &N, &R, &R, &one, exxiyV, &N, const_cast<double *> (Dy), &R, &zero, Vtemp, &N);
	dgemm_(transn, transn, &N, &R, &R, &one, Vtemp, &N, const_cast<double *> (Dy), &R, &zero, exxiyV, &N);
	dgemm_(transn, transt, &N, &R, &R, &one, const_cast<double *> (Vy), &N, exxiyD, &R, &zero, Vtemp, &N);
	dgemm_(transn, transn, &N, &R, &R, &one, Vtemp, &N, const_cast<double *> (Dy), &R, &zero, VxiDtDy, &N);
	daxpy_(&NR, &one, VxiDtDy, &inc, exxiyV, &inc);

	delete[] VxiDtDy;

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	ProdVector *prodexetax = dynamic_cast<ProdVector *> (exetax);
	ProdVector *prodexxiy = dynamic_cast<ProdVector *> (exxiy);
	ProdVector *prodexresult = dynamic_cast<ProdVector *> (exresult);
	manifolds[0]->ExtrProjection(prody->GetElement(0), prodexxiy->GetElement(0), prodexxiy->GetElement(0));
	manifolds[2]->ExtrProjection(prody->GetElement(2), prodexxiy->GetElement(2), prodexxiy->GetElement(2));
	prodexresult->NewMemoryOnWrite();

	ProductManifold::coTangentVector(x, exetax, y, exxiy, exresult);

	ExtrProjectionStiePerp(prodx->GetElement(0), prodexresult->GetElement(0), prodexresult->GetElement(0));
	ExtrProjectionStiePerp(prodx->GetElement(2), prodexresult->GetElement(2), prodexresult->GetElement(2));


	//- It is not a good approach -------
	LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
	const double *Dx = LRx->GetElement(1)->ObtainReadData();
	double *dU = prodexresult->GetElement(0)->ObtainWritePartialData();
	double *dV = prodexresult->GetElement(2)->ObtainWritePartialData();
	integer RR = R * R, info;
	integer *IPIV = new integer[r];
	double *Dinv = new double[r * r + m * r + n * r];
	Utemp = Dinv + r * r;
	Vtemp = Utemp + m * r;
	dcopy_(&RR, const_cast<double *> (Dx), &inc, Dinv, &inc);
	dgetrf_(&R, &R, Dinv, &R, IPIV, &info);
	integer lwork = -1;
	double lworkopt;
	dgetri_(&R, Dinv, &R, IPIV, &lworkopt, &lwork, &info);
	lwork = static_cast<integer> (lworkopt);
	double *work = new double[lwork];
	dgetri_(&R, Dinv, &R, IPIV, work, &lwork, &info);
	delete[] work;
	delete[] IPIV;

	dgemm_(transn, transt, &M, &R, &R, &one, dU, &M, Dinv, &R, &zero, Utemp, &M);
	dgemm_(transn, transn, &M, &R, &R, &one, Utemp, &M, Dinv, &R, &zero, dU, &M);

	dgemm_(transn, transn, &N, &R, &R, &one, dV, &N, Dinv, &R, &zero, Vtemp, &N);
	dgemm_(transn, transt, &N, &R, &R, &one, Vtemp, &N, Dinv, &R, &zero, dV, &N);

	delete[] Dinv;
	ObtainIntr(x, exresult, result);
	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i]->SetIsIntrApproach(true);
	}
	delete exetax;
	delete exxiy;
	delete exresult;
};

void LowRank::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
{
	Vector *exetax = EMPTYEXTR->ConstructEmpty();
	Vector *exxix = EMPTYEXTR->ConstructEmpty();
	Vector *exresult = EMPTYEXTR->ConstructEmpty();
	ObtainExtr(x, etax, exetax);
	ObtainExtr(x, xix, exxix);

	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i]->SetIsIntrApproach(false);
	}

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodexetax = dynamic_cast<ProdVector *> (exetax);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	ProdVector *prodexxix = dynamic_cast<ProdVector *> (exxix);
	ProdVector *prodexresult = dynamic_cast<ProdVector *> (exresult);
	prodexresult->NewMemoryOnWrite();
	manifolds[0]->DiffRetraction(prodx->GetElement(0), prodexetax->GetElement(0), prody->GetElement(0),
		prodexxix->GetElement(0), prodexresult->GetElement(0), IsEtaXiSameDir);
	manifolds[1]->DiffRetraction(prodx->GetElement(1), prodexetax->GetElement(1), prody->GetElement(1),
		prodexxix->GetElement(1), prodexresult->GetElement(1), IsEtaXiSameDir);
	manifolds[2]->DiffRetraction(prodx->GetElement(2), prodexetax->GetElement(2), prody->GetElement(2),
		prodexxix->GetElement(2), prodexresult->GetElement(2), IsEtaXiSameDir);

	ObtainIntr(y, exresult, result);
	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i]->SetIsIntrApproach(true);
	}
	delete exetax;
	delete exxix;
	delete exresult;

	if (IsEtaXiSameDir)
	{
		const double *etaxTV = etax->ObtainReadData();
		const double *xixTV = xix->ObtainReadData();
		double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
		SharedSpace *beta = new SharedSpace(1, 1);
		double *betav = beta->ObtainWriteEntireData();
		betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
		etax->AddToTempData("beta", beta);

		Vector *TReta = result->ConstructEmpty();
		result->CopyTo(TReta);
		ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
		SharedSpace *SharedTReta = new SharedSpace(TReta);
		etax->AddToTempData("betaTReta", SharedTReta);
	}
};

void LowRank::ExtrProjection(Variable *x, Vector *etax, Vector *result) const
{
	Vector *inetax = EMPTYINTR->ConstructEmpty();
	ObtainIntr(x, etax, inetax);
	ObtainExtr(x, inetax, result);
	delete inetax;
};

void LowRank::ExtrProjectionStiePerp(Variable *x, Vector *v, Vector *result) const
{
	integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
	double *UtV = new double[P * P];
	const double *U = x->ObtainReadData();
	const double *V = v->ObtainReadData();
	double *resultTV = result->ObtainWriteEntireData();

	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (U), &N, const_cast<double *> (V), &N, &zero, UtV, &P);

	if (V != resultTV)
		dcopy_(&Length, const_cast<double *> (V), &inc, resultTV, &inc);
	double negone = -1;
	dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (U), &N, UtV, &P, &one, resultTV, &N);
	delete[] UtV;
};

void LowRank::SetParams(PARAMSMAP params)
{
	Manifold::SetParams(params);
};
