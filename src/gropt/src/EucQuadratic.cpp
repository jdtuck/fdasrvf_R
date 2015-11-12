
#include "EucQuadratic.h"

EucQuadratic::EucQuadratic(double *M, integer dim)
{
	Dim = dim;
	A = M;
};

EucQuadratic::~EucQuadratic(void)
{
};

double EucQuadratic::f(Variable *x) const
{
	const double *v = x->ObtainReadData();
	SharedSpace *Temp = new SharedSpace(1, Dim);
	double *temp = Temp->ObtainWriteEntireData();

	char *transn = const_cast<char *> ("n");
	double one = 1, zero = 0;
	integer inc = 1, N = Dim;
	dgemv_(transn, &N, &N, &one, A, &N, const_cast<double *> (v), &inc, &zero, temp, &inc);

	x->AddToTempData("Ax", Temp);

	return ddot_(&N, const_cast<double *> (v), &inc, temp, &inc);
};

void EucQuadratic::Grad(Variable *x, Vector *gf) const
{
	double *gfTV = gf->ObtainWriteEntireData();
	const SharedSpace *Temp = x->ObtainReadTempData("Ax");
	const double *v = Temp->ObtainReadData();

	integer N = Dim, inc = 1;
	dcopy_(&N, const_cast<double *> (v), &inc, gfTV, &inc);
	double two = 2;
	dscal_(&N, &two, gfTV, &inc);
};

void EucQuadratic::HessianEta(Variable *x, Vector *etax, Vector *xix) const
{
	const double *v = etax->ObtainReadData();
	double *xixTV = xix->ObtainWriteEntireData();

	char *transn = const_cast<char *> ("n");
	integer N = Dim, inc = 1;
	double two = 2, zero = 0;
	dgemv_(transn, &N, &N, &two, A, &N, const_cast<double *> (v), &inc, &zero, xixTV, &inc);
};
