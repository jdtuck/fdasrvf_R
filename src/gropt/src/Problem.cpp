
#include "Problem.h"

void Problem::CheckGradHessian(const Variable *xin) const
{
	UseGrad = true;
	UseHess = true;
	integer length;
	double normxi;
	double t, fx, fy;
	double *X, *Y;
	Vector *etax;
	Variable *x = xin->ConstructEmpty();
	xin->CopyTo(x);
	if (Domain->GetIsIntrinsic())
		etax = Domain->GetEMPTYINTR()->ConstructEmpty();
	else
		etax = Domain->GetEMPTYEXTR()->ConstructEmpty();
	etax->RandUnform();
	Vector *xi = etax->ConstructEmpty();
	Vector *gfx = etax->ConstructEmpty();
	Vector *Hv = etax->ConstructEmpty();
	Variable *y = x->ConstructEmpty();
	fx = f(x);
	Grad(x, gfx);
	gfx->CopyTo(etax);//--
	//double *etaxTV = etax->ObtainWriteEntireData();///---
	//integer nnn = etax->Getlength();
	//for (integer i = 0; i < nnn; i++)//--
	//{
	//	etaxTV[i] = sin(static_cast<double> (i) / (etax->Getlength() - 1) / 2);
	//}
	//for (integer i = 0; i < 5; i++)//---
	//	etaxTV[nnn - 1 - i] = 0;//--

	//etax->Print("etax:");//--
	Domain->Projection(x, etax, xi);
	normxi = sqrt(Domain->Metric(x, xi, xi));
	Domain->ScaleTimesVector(x, 100.0 / normxi, xi, xi); // initial length of xi is 100
	//xi->Print("xi:");//---
	// the length of xi variances from 100 to 100*2^(-35) approx 6e-9
	t = 1;
	length = 35;
	X = new double [length * 2];
	Y = X + length;
	for (integer i = 0; i < length; i++)
	{
		Domain->Retraction(x, xi, y);
		fy = f(y);
		HessianEta(x, xi, Hv);
		Y[i] = log(fabs(fy - fx - Domain->Metric(x, gfx, xi) - 0.5 * Domain->Metric(x, xi, Hv)));
		X[i] = 0.5 * log(Domain->Metric(x, xi, xi));
		printf("i:%d,|eta|:%.3e,(fy-fx)/<gfx,eta>:%.3e,(fy-fx-<gfx,eta>)/<0.5 eta, Hessian eta>:%.3e\n", i, 
			sqrt(Domain->Metric(x, xi, xi)), (fy-fx)/Domain->Metric(x, gfx, xi), 
			(fy - fx - Domain->Metric(x, gfx, xi)) / (0.5 * Domain->Metric(x, xi, Hv)));
		Domain->ScaleTimesVector(x, 0.5, xi, xi);
	}

	std::cout << "CHECK GRADIENT:" << std::endl;
	std::cout << "\tSuppose the point is not a critical point." << std::endl;
	std::cout << "\tIf there exists an interval of |eta| such that (fy - fx) / <gfx, eta>" << std::endl;
	std::cout << "\tapproximates ONE, then the gradient is probably correct!" << std::endl;

	std::cout << "CHECK THE ACTION OF THE HESSIAN (PRESUME GRADIENT IS CORRECT):" << std::endl;
	std::cout << "\tSuppose the retraction is second order or the point is a critical point." << std::endl;
	std::cout << "\tIf there exists an interval of |eta| such that (fy-fx-<gfx,eta>)/<0.5 eta, Hessian eta>" << std::endl;
	std::cout << "\tapproximates ONE, then the action of Hessian is probably correct." << std::endl;

	////TEST IDEA2: 
	//for (integer i = 1; i < length - 1; i++)
	//	printf("log(|eta|):%.3e, slope:%.3e\n", X[i], (Y[i + 1] - Y[i - 1]) / (X[i + 1] - X[i - 1]));
	//std::cout << "CHECK GRADIENT:" << std::endl;
	//std::cout << "\tIf there exists an interval of |eta| such that the slopes " << std::endl;
	//std::cout << "\tapproximate TWO, then the gradient is probably correct!" << std::endl;

	//std::cout << "CHECK THE ACTION OF THE HESSIAN (PRESUME GRADIENT IS CORRECT AND" << std::endl;
	//std::cout << "THE COST FUNCTION IS NOT ONLY QUADRATIC):" << std::endl;
	//std::cout << "\tIf there exists an interval of |eta| such that the slopes" << std::endl;
	//std::cout << "\tapproximate THREE, then the action of Hessian is probably correct." << std::endl;

	//x->Print("1, x:", false);//---
	delete xi;
	//x->Print("2, x:", false);//---
	//gfx->Print("2, gfx:", false);//---
	delete gfx;
	//x->Print("3, x:", false);//---
	delete y;
	//x->Print("4, x:", false);//---
	delete Hv;
	//x->Print("5, x:", false);//---
	delete[] X;
	delete etax;
	//x->Print("x:", false);//---
	delete x;
};

void Problem::Grad(Variable *x, Vector *gf) const
{
	if (!Domain->GetIsIntrinsic())
	{
		RieGrad(x, gf);
		return;
	}
	Vector *exgf = Domain->GetEMPTYEXTR()->ConstructEmpty();
	RieGrad(x, exgf);
	//exgf->Print("exgf:");//---
	Domain->ObtainIntr(x, exgf, gf);
	delete exgf;
};

void Problem::HessianEta(Variable *x, Vector *etax, Vector *xix) const
{
	if (!Domain->GetIsIntrinsic())
	{
		RieHessianEta(x, etax, xix);
		return;
	}

	Vector *exxix = Domain->GetEMPTYEXTR()->ConstructEmpty();
	Vector *exetax = Domain->GetEMPTYEXTR()->ConstructEmpty();
	Domain->ObtainExtr(x, etax, exetax);
	RieHessianEta(x, exetax, exxix);
	Domain->ObtainIntr(x, exxix, xix);
	delete exxix;
	delete exetax;
};

void Problem::RieGrad(Variable *x, Vector *gf) const
{
	EucGrad(x, gf);
	Domain->EucGradToGrad(x, gf, gf, this);
};

void Problem::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const 
{
	EucHessianEta(x, etax, xix);
	Domain->EucHvToHv(x, etax, xix, xix, this);
};

void Problem::EucGrad(Variable *x, Vector *egf) const
{
	std::cout << "Euclidean Gradient has not been done!" << std::endl;
};

void Problem::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
{
	std::cout << "The action of Euclidean Hessian has not been done!" << std::endl;
};

void Problem::SetDomain(Manifold *inDomain)
{
	Domain = inDomain;
};

Problem::~Problem(void)
{
};

void Problem::SetUseGrad(bool usegrad) const
{
	UseGrad = usegrad;
};

void Problem::SetUseHess(bool usehess) const
{
	UseHess = usehess;
};
