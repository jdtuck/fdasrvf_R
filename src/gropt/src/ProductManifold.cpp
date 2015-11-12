#include "ProductManifold.h"

ProductManifold::ProductManifold(integer numberofmanifolds, ...)
{
	numofmani = numberofmanifolds;
	powsinterval = new integer[numofmani + 1];
	manifolds = new Manifold *[numofmani];
	//integer numofargu = 2 * numberofmanifolds;
	va_list argptr;
	va_start(argptr, numberofmanifolds);
	powsinterval[0] = 0;
	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i] = va_arg(argptr, Manifold *);
		powsinterval[i + 1] = powsinterval[i] + va_arg(argptr, integer);
	}
	va_end(argptr);

	HasHHR = false;
	HasLockCon = false;
	numoftotalmani = 0;
	ExtrinsicDim = 0;
	IntrinsicDim = 0;
	for (integer i = 0; i < numofmani; i++)
	{
		ExtrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetExtrDim();
		IntrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetIntrDim();
		numoftotalmani += (powsinterval[i + 1] - powsinterval[i]);
	}
	name.assign("Product Manifold");
	// For product manifold, individual manifolds are using their parameters, intrapproach, 
	// to specify whether intrinsic approach is used or not.
	IsIntrApproach = true;
	//for (integer i = 0; i < numofmani; i++)
	//{
	//	manifolds[i]->SetIsIntrApproach(true);
	//}

	Element **elements = new Element *[numoftotalmani];
	for (integer i = 0; i < numofmani; i++)
	{
		if (manifolds[i]->GetIsIntrinsic())
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYINTR());
			}
		}
		else
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYEXTR());
			}
		}
	}
	// EMPTYINTR specify the tangent vector using the approach given by the IsIntrinsic parameters of each manfold.
	EMPTYINTR = new ProductElement(elements, numoftotalmani, powsinterval, numofmani);

	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYEXTR());
		}
	}
	// EMPTYEXTR use extrinsic approach for all the manifolds
	EMPTYEXTR = new ProductElement(elements, numoftotalmani, powsinterval, numofmani);
	delete[] elements;
};

ProductManifold::ProductManifold(Manifold **inmanifolds, integer innumofmani, integer *inpowsinterval, integer innumoftotalmani)
{
	numofmani = innumofmani;
	powsinterval = new integer[numofmani + 1];
	manifolds = new Manifold *[numofmani];
	powsinterval[0] = 0;
	for (integer i = 0; i < numofmani; i++)
	{
		manifolds[i] = inmanifolds[i];
		powsinterval[i + 1] = inpowsinterval[i + 1];
	}

	HasHHR = false;
	HasLockCon = false;
	numoftotalmani = 0;
	ExtrinsicDim = 0;
	IntrinsicDim = 0;
	for (integer i = 0; i < numofmani; i++)
	{
		ExtrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetExtrDim();
		IntrinsicDim += (powsinterval[i + 1] - powsinterval[i]) * manifolds[i]->GetIntrDim();
		numoftotalmani += (powsinterval[i + 1] - powsinterval[i]);
	}
	name.assign("Product Manifold");
	// For product manifold, either all the individual manifolds are using intrinsic approach
	// or all of them are using extrinsic approach.
	IsIntrApproach = true;
	//for (integer i = 0; i < numofmani; i++)
	//{
	//	manifolds[i]->SetIsIntrApproach(true);
	//}

	Element **elements = new Element *[numoftotalmani];
	for (integer i = 0; i < numofmani; i++)
	{
		if (manifolds[i]->GetIsIntrinsic())
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYINTR());
			}
		}
		else
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYEXTR());
			}
		}
	}
	// EMPTYINTR specify the tangent vector using the approach given by the IsIntrinsic parameters of each manfold.
	EMPTYINTR = new ProductElement(elements, numoftotalmani, powsinterval, numofmani);

	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYEXTR());
		}
	}
	// EMPTYEXTR use extrinsic approach for all the manifolds
	EMPTYEXTR = new ProductElement(elements, numoftotalmani, powsinterval, numofmani);
	delete[] elements;
/*
	Element **elements = new Element *[numoftotalmani];
	for (integer i = 0; i < numofmani; i++)
	{
		if (manifolds[i]->GetIsIntrinsic())
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYINTR());
			}
		}
		else
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYEXTR());
			}
		}
	}
	EMPTYINTR = new ProductElement(elements, numoftotalmani, powsinterval, numofmani);
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			elements[j] = const_cast<Element *> (manifolds[i]->GetEMPTYEXTR());
		}
	}
	EMPTYEXTR = new ProductElement(elements, numoftotalmani, powsinterval, numofmani);
	delete[] elements;*/
};

ProductManifold::~ProductManifold(void)
{
	delete EMPTYINTR;
	delete EMPTYEXTR;
	delete[] manifolds;
	delete[] powsinterval;
};

double ProductManifold::Metric(Variable *x, Vector *etax, Vector *xix) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVector *prodxix = dynamic_cast<ProdVector *> (xix);
	double result = 0;
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			result += manifolds[i]->Metric(prodx->GetElement(j), prodetax->GetElement(j), prodxix->GetElement(j));
		}
	}

	return result;
};

void ProductManifold::LinearOPEEta(Variable *x, LinearOPE *Hx, Vector *etax, Vector *result) const
{
	Manifold::LinearOPEEta(x, Hx, etax, result);

#ifdef CHECKMEMORY
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	prodresult->CheckMemory();
#endif
};

void ProductManifold::ScaleTimesVector(Variable *x, double scaler, Vector *etax, Vector *result) const
{
	Manifold::ScaleTimesVector(x, scaler, etax, result);

#ifdef CHECKMEMORY
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	prodresult->CheckMemory();
#endif
};

void ProductManifold::VectorLinearCombination(Variable *x, double scaler1, Vector *etax, double scaler2, Vector *xix, Vector *result) const
{
	Manifold::VectorLinearCombination(x, scaler1, etax, scaler2, xix, result);

#ifdef CHECKMEMORY
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	prodresult->CheckMemory();
#endif
};

void ProductManifold::Projection(Variable *x, Vector *v, Vector *result) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodv = dynamic_cast<ProdVector *> (v);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	if (v == result)
	{
		ProdVector *prodresultTemp = prodresult->ConstructEmpty();
		prodresultTemp->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->Projection(prodx->GetElement(j), prodv->GetElement(j), prodresultTemp->GetElement(j));
			}
		}
		prodresultTemp->CopyTo(result);
		delete prodresultTemp;
	}
	else
	{
		prodresult->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->Projection(prodx->GetElement(j), prodv->GetElement(j), prodresult->GetElement(j));
			}
		}
	}

#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::RandomTangentVectors(Variable *x, integer N, Vector **result_arr) const
{
};

void ProductManifold::Retraction(Variable *x, Vector *etax, Variable *result) const
{
#ifdef TESTELASTICCURVESRO
	if (x->TempDataExist("w"))
	{
		const SharedSpace *Sharedw = x->ObtainReadTempData("w");
		const double *wptr = Sharedw->ObtainReadData();
		SharedSpace *Sharedww = new SharedSpace(1, 1);
		double *wwptr = Sharedww->ObtainWriteEntireData();
		wwptr[0] = wptr[0] * 0.8;
		result->AddToTempData("w", Sharedww);
	}
#endif

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	prodresult->NewMemoryOnWrite();
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			manifolds[i]->Retraction(prodx->GetElement(j), prodetax->GetElement(j), prodresult->GetElement(j));
		}
	}

#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	ProdVector *prodxiy = dynamic_cast<ProdVector *> (xiy);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);

	if (xiy == result)
	{
		ProdVector *prodresultTemp = prodresult->ConstructEmpty();
		prodresultTemp->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->coTangentVector(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxiy->GetElement(j), prodresultTemp->GetElement(j));
			}
		}
		prodresultTemp->CopyTo(prodresult);
		delete prodresultTemp;
	}
	else
	{
		prodresult->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->coTangentVector(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxiy->GetElement(j), prodresult->GetElement(j));
			}
		}
	}
#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

double ProductManifold::Beta(Variable *x, Vector *etax) const
{
	if (!HasHHR)
		return 1;

	if (etax->TempDataExist("beta"))
	{
		const SharedSpace *beta = etax->ObtainReadTempData("beta");
		const double *betav = beta->ObtainReadData();
		return betav[0];
	}

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	Variable *z;
	const SharedSpace *beta;
	const double *betav;
	double numerator = 0, denominator = 0;
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			z = prodx->GetElement(j);
			if (z->TempDataExist("beta"))
			{
				beta = z->ObtainReadTempData("beta");
				betav = beta->ObtainReadData();
				numerator += betav[1];
				denominator += betav[2];
			}
			else
			{
				numerator += manifolds[j]->Metric(prodx->GetElement(j), prodetax->GetElement(j), prodetax->GetElement(j));
				denominator += numerator;
			}
		}
	}
	return sqrt(numerator / denominator);
};

void ProductManifold::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	ProdVector *prodxix = dynamic_cast<ProdVector *> (xix);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	if (xix == result)
	{
		ProdVector *prodresultTemp = prodresult->ConstructEmpty();
		prodresultTemp->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->DiffRetraction(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxix->GetElement(j), prodresultTemp->GetElement(j), IsEtaXiSameDir);
			}
		}
		prodresultTemp->CopyTo(prodresult);
		delete prodresultTemp;
	}
	else
	{
		prodresult->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->DiffRetraction(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxix->GetElement(j), prodresult->GetElement(j), IsEtaXiSameDir);
			}
		}
	}

#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif

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

void ProductManifold::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
{
	if (HasHHR)
		return LCVectorTransport(x, etax, y, xix, result);

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	ProdVector *prodxix = dynamic_cast<ProdVector *> (xix);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);

	if (xix == result)
	{
		ProdVector *prodresultTemp = prodresult->ConstructEmpty();
		prodresultTemp->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->VectorTransport(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxix->GetElement(j), prodresultTemp->GetElement(j));
			}
		}
		prodresultTemp->CopyTo(prodresult);
		delete prodresultTemp;
	}
	else
	{
		prodresult->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->VectorTransport(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxix->GetElement(j), prodresult->GetElement(j));
			}
		}
	}
#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
{
	if (HasHHR)
		return LCInverseVectorTransport(x, etax, y, xiy, result);

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	ProdVector *prodxiy = dynamic_cast<ProdVector *> (xiy);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);

	if (xiy == result)
	{
		ProdVector *prodresultTemp = prodresult->ConstructEmpty();
		prodresultTemp->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->InverseVectorTransport(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxiy->GetElement(j), prodresultTemp->GetElement(j));
			}
		}
		prodresultTemp->CopyTo(prodresult);
		delete prodresultTemp;
	}
	else
	{
		prodresult->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->InverseVectorTransport(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), prodxiy->GetElement(j), prodresult->GetElement(j));
			}
		}
	}
#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
{
	if (HasHHR)
		return LCTranHInvTran(x, etax, y, Hx, result);

	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
	integer start, end = 0;
	Hx->CopyTo(result);
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			start = end;
			end = start + prodetax->GetElement(j)->Getlength();
			manifolds[i]->HInvTran(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), result, start, end, result);
			manifolds[i]->TranH(prodx->GetElement(j), prodetax->GetElement(j), prody->GetElement(j), result, start, end, result);
		}
	}
};

void ProductManifold::HaddScaledRank1OPE(Variable *x, LinearOPE *Hx, double scaler, Vector *etax, Vector *xix, LinearOPE *result) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodxix = dynamic_cast<ProdVector *> (xix);

	ProdVector *prodxixflat = prodxix->ConstructEmpty();
	prodxixflat->NewMemoryOnWrite();

	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			manifolds[i]->ObtainEtaxFlat(prodx->GetElement(j), prodxix->GetElement(j), prodxixflat->GetElement(j));
		}
	}
	Manifold::HaddScaledRank1OPE(x, Hx, scaler, etax, prodxixflat, result);
	delete prodxixflat;
};

void ProductManifold::ProdElementToElement(const ProductElement *ProdElem, Element *Elem) const
{
	const double *Prodspace;
	double *space = Elem->ObtainWriteEntireData();
	integer N, inc = 1, idx = 0;
	for (integer i = 0; i < numoftotalmani; i++)
	{
		Prodspace = ProdElem->GetElement(i)->ObtainReadData();
		N = ProdElem->GetElement(i)->Getlength();
		dcopy_(&N, const_cast<double *> (Prodspace), &inc, space + idx, &inc);
		idx += N;
	}
};

void ProductManifold::ElementToProdElement(const Element *Elem, ProductElement *ProdElem) const
{
	double *Prodspace;
	const double *space = Elem->ObtainReadData();
	integer N, inc = 1, idx = 0;
	for (integer i = 0; i < numoftotalmani; i++)
	{
		Prodspace = ProdElem->GetElement(i)->ObtainWriteEntireData();
		N = ProdElem->GetElement(i)->Getlength();
		dcopy_(&N, const_cast<double *> (space + idx), &inc, Prodspace, &inc);
		idx += N;
	}
};

void ProductManifold::ExtrProjection(Variable *x, Vector *v, Vector *result) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodv = dynamic_cast<ProdVector *> (v);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	if (v == result)
	{
		ProdVector *prodresultTemp = prodresult->ConstructEmpty();
		prodresultTemp->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->ExtrProjection(prodx->GetElement(j), prodv->GetElement(j), prodresultTemp->GetElement(j));
			}
		}
		prodresultTemp->CopyTo(prodresult);
		delete prodresultTemp;
	}
	else
	{
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->ExtrProjection(prodx->GetElement(j), prodv->GetElement(j), prodresult->GetElement(j));
			}
		}
	}
#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	prodresult->NewMemoryOnWrite();
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			if (manifolds[i]->GetIsIntrinsic())
			{
				manifolds[i]->ObtainIntr(prodx->GetElement(j), prodetax->GetElement(j), prodresult->GetElement(j));
			}
			else
			{
				prodetax->GetElement(j)->CopyTo(prodresult->GetElement(j));
			}
		}
	}
#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodintretax = dynamic_cast<ProdVector *> (intretax);
	ProdVector *prodresult = dynamic_cast<ProdVector *> (result);
	prodresult->NewMemoryOnWrite();
	for (integer i = 0; i < numofmani; i++)
	{
		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			if (manifolds[i]->GetIsIntrinsic())
			{
				manifolds[i]->ObtainExtr(prodx->GetElement(j), prodintretax->GetElement(j), prodresult->GetElement(j));
			}
			else
			{
				prodintretax->GetElement(j)->CopyTo(prodresult->GetElement(j));
			}
		}
	}
#ifdef CHECKMEMORY
	prodresult->CheckMemory();
#endif
};

void ProductManifold::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodegf = dynamic_cast<ProdVector *> (egf);
	ProdVector *prodgf = dynamic_cast<ProdVector *> (gf);
	if (egf == gf)
	{
		ProdVector *prodgfTemp = prodgf->ConstructEmpty();
		prodgfTemp->NewMemoryOnWrite();

		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->EucGradToGrad(prodx->GetElement(j), prodegf->GetElement(j), prodgfTemp->GetElement(j), prob);
			}
		}
		prodgfTemp->CopyTo(prodgf);
		delete prodgfTemp;
	}
	else
	{
		prodgf->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->EucGradToGrad(prodx->GetElement(j), prodegf->GetElement(j), prodgf->GetElement(j), prob);
			}
		}
	}
#ifdef CHECKMEMORY
	prodgf->CheckMemory();
#endif
};

void ProductManifold::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector *xix, const Problem *prob) const
{
	ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
	ProdVector *prodetax = dynamic_cast<ProdVector *> (etax);
	ProdVector *prodexix = dynamic_cast<ProdVector *> (exix);
	ProdVector *prodxix = dynamic_cast<ProdVector *> (xix);
	if (exix == xix)
	{
		ProdVector *prodxixTemp = prodxix->ConstructEmpty();
		prodxixTemp->NewMemoryOnWrite();

		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->EucHvToHv(prodx->GetElement(j), prodetax->GetElement(j), prodexix->GetElement(j), prodxixTemp->GetElement(j), prob);
			}
		}
		prodxixTemp->CopyTo(prodxix);
		delete prodxixTemp;
	}
	else
	{
		prodxix->NewMemoryOnWrite();
		for (integer i = 0; i < numofmani; i++)
		{
			for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
			{
				manifolds[i]->EucHvToHv(prodx->GetElement(j), prodetax->GetElement(j), prodexix->GetElement(j), prodxix->GetElement(j), prob);
			}
		}
	}
#ifdef CHECKMEMORY
	prodxix->CheckMemory();
#endif
};

void ProductManifold::CheckParams(void) const
{
	if (numoftotalmani == 1)
	{
		manifolds[0]->CheckParams();
	}
	else
	{
		Manifold::CheckParams();
		for (integer i = 0; i < numofmani; i++)
		{
			std::cout << i << "-th manifold parameters (the number is " << powsinterval[i + 1] - powsinterval[i] << ") :" << std::endl;
			manifolds[i]->CheckParams();
		}
	}
};
