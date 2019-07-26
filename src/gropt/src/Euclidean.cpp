
#include "Euclidean.h"

/*Define the namespace*/
namespace ROPTLIB{

	Euclidean::Euclidean(integer r, integer c, integer n)
	{
		row = r;
		col = c;
		num = n;
		IsIntrApproach = false;
		HasHHR = false;
		UpdBetaAlone = false;
		name.assign("Euclidean");
		IntrinsicDim = n * r * c;
		ExtrinsicDim = n * r * c;
		EMPTYEXTR = new EucVector(r, c, n);
		EMPTYINTR = new EucVector(r, c, n);
	};

	Euclidean::~Euclidean(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	void Euclidean::CheckParams(void) const
	{
		Manifold::CheckParams();
		Rprintf("%s PARAMETERS:\n", name.c_str());
		if (col == 1 && num == 1)
			Rprintf("row           :%15d\n", row);
		else
			if (num == 1)
			{
				Rprintf("row           :%15d,\t", row);
				Rprintf("col           :%15d\n", col);
			}
			else
			{
				Rprintf("row           :%15d,\t", row);
				Rprintf("col           :%15d\n", col);
				Rprintf("num           :%15d\n", num);
			}
	};

	void Euclidean::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		egf->CopyTo(gf);
	};

	void Euclidean::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		exix->CopyTo(xix);
	};
}; /*end of ROPTLIB namespace*/
