
#ifndef ELASTICCURVESRO_H
#define ELASTICCURVESRO_H

//#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)

#include <L2Sphere.h>
#include <L2SphereVariable.h>
#include <L2SphereVector.h>
#include <Euclidean.h>
#include <EucVariable.h>
#include <EucVector.h>
#include <OrthGroup.h>
#include <OrthGroupVariable.h>
#include <OrthGroupVector.h>
#include <ProductElement.h>
#include <ProductManifold.h>

#include <def.h>
#include <Problem.h>
#include <SharedSpace.h>
#include <Spline.h>

#define NNBRS 23
const int Nbrs[NNBRS][2] = {
	{ 1, 1 },
	{ 1, 2 },
	{ 2, 1 },
	{ 2, 3 },
	{ 3, 2 },
	{ 1, 3 },
	{ 3, 1 },
	{ 1, 4 },
	{ 3, 4 },
	{ 4, 3 },
	{ 4, 1 },
	{ 1, 5 },
	{ 2, 5 },
	{ 3, 5 },
	{ 4, 5 },
	{ 5, 4 },
	{ 5, 3 },
	{ 5, 2 },
	{ 5, 1 },
	{ 1, 6 },
	{ 5, 6 },
	{ 6, 5 },
	{ 6, 1 }
};

// Find best reparameterization (rotation) for functions, open and closed curves in R^d
// \int_D \|Oq1(t) - q2(\int_0^t l^4(s) ds + m mod 1) l^2(t)\|_2^2 dt
// the domain D = [0, 1] for open curves, D is a circle for closed curves.
class ElasticCurvesRO : public Problem{
public:
	// if isclosed is true, then make sure the first and last points of C1 and C2 must be the same.
	ElasticCurvesRO(double *inq1, double *inq2, integer ind, integer inn, double inw, bool inrotated, bool inisclosed);
	virtual ~ElasticCurvesRO();

	//virtual void SetDomain(Manifold *inDomain);
//	Element **GetInitialXandSetDomain(integer &length);

	virtual double f(Variable *x) const;

	virtual void EucGrad(Variable *x, Vector *egf) const;
	virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	static void PointwiseInnerProd(const double *q1, const double *q2, integer d, integer n, double *result);
	static void PointwiseQProdl(const double *q1, const double *l, integer d, integer n, double *result);
	static void PointwiseProd(const double *l1, const double *l2, integer n, double *result);
	static void CumTrapz(const double *l, integer n, double intv, double *result);
	static double Trapz(const double *l, integer n, double intv);

	mutable double w;			// coefficient of penlty term
private:

	double *q1;			// n by d matrices
	double *q2_coefs;	// (n - 1) by 4 d matrices, d cubic splines
	double *dq2_coefs;	// (n - 1) by 3 d matrices, d cubic splines
	double *ddq2_coefs;	// (n - 1) by 2 d matrices, d cubic splines
	mutable integer n;			// the number of points for representing a function/curve
	mutable integer d;			// the dimension of the curves are
	bool rotated;		// involve rotation or not
	bool isclosed;		// for closed curves or open curves
};

#endif // ELASTICCURVESRO_H
