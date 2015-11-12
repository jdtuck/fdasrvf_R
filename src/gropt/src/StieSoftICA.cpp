
#include "StieSoftICA.h"

StieSoftICA::StieSoftICA(double *inCs, integer inn, integer inp, integer inN)
{
	Cs = inCs;
	n = inn;
	p = inp;
	N = inN;
};

double StieSoftICA::f(Variable *x) const
{
	const double *xxM = x->ObtainReadData();
	SharedSpace *SharedCY = new SharedSpace(1, n * p * N);
	SharedSpace *SharedD = new SharedSpace(1, p * N);
	double *CY = SharedCY->ObtainWriteEntireData();
	double *D = SharedD->ObtainWriteEntireData();

	char *transn = const_cast<char *> ("n");
	double one = 1, zero = 0;
	integer inc = 1;
	for (integer i = 0; i < N; i++)
	{
		dgemm_(transn, transn, &n, &p, &n, &one, Cs + n * n * i, &n, const_cast<double *> (xxM), &n, &zero, CY + n * p * i, &n);
	}

	for (integer i = 0; i < N; i++)
	{
		for (integer j = 0; j < p; j++)
		{
			D[j + i * p] = ddot_(&n, const_cast<double *> (xxM + j * n), &inc, CY + n * p * i + j * n, &inc);
		}
	}
	integer pN = p * N;
	double result = - ddot_(&pN, D, &inc, D, &inc);
	if (UseGrad)
	{
		x->AddToTempData("CY", SharedCY);
		x->AddToTempData("D", SharedD);
	}
	else
	{
		delete SharedCY;
		delete SharedD;
	}
	return result;
};

void StieSoftICA::EucGrad(Variable *x, Vector *egf) const
{
	const SharedSpace *SharedCY = x->ObtainReadTempData("CY");
	const double *CY = SharedCY->ObtainReadData();
	const SharedSpace *SharedD = x->ObtainReadTempData("D");
	const double *D = SharedD->ObtainReadData();

	double *egfTV = egf->ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
	{
		egfTV[i] = 0;
	}
	integer np = n * p, inc = 1;
	double coef = 0;
	for (integer i = 0; i < N; i++)
	{
		for (integer j = 0; j < p; j++)
		{
			coef = - D[j + i * p] * 4;
			daxpy_(&n, &coef, const_cast<double *> (CY + i * n * p + j * n), &inc, egfTV + j * n, &inc);
		}
	}
};

void StieSoftICA::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
{
	const double *xxM = x->ObtainReadData();
	const double *etaxTV = etax->ObtainReadData();
	double *exixTV = exix->ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
	{
		exixTV[i] = 0;
	}

	const SharedSpace *SharedCY = x->ObtainReadTempData("CY");
	const double *CY = SharedCY->ObtainReadData();
	const SharedSpace *SharedD = x->ObtainReadTempData("D");
	const double *D = SharedD->ObtainReadData();

	double *temp = new double[n * p];
	double coef = 0;
	integer np = n * p, inc = 1;
	char *transn = const_cast<char *> ("n");
	double one = 1, zero = 0;
	for (integer i = 0; i < N; i++)
	{
		dcopy_(&np, const_cast<double *> (etaxTV), &inc, temp, &inc);
		for (integer j = 0; j < p; j++)
		{
			dscal_(&n, const_cast<double *> (D + i * p + j), temp + j * n, &inc);
		}
		for (integer j = 0; j < p; j++)
		{
			coef = ddot_(&n, const_cast<double *> (etaxTV + j * n), &inc, const_cast<double *> (CY + i * n * p + j * n), &inc) * 2;
			daxpy_(&n, &coef, const_cast<double *> (xxM + j * n), &inc, temp + j * n, &inc);
		}

		dgemm_(transn, transn, &n, &p, &n, &one, Cs + i * n * n, &n, temp, &n, &one, exixTV, &n);
	}

	delete[] temp;

	Domain->ScaleTimesVector(x, -4.0, exix, exix);
};
