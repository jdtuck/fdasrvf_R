
#include "Euclidean.h"

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
	std::cout << name << " PARAMETERS:" << std::endl;
	if (col == 1 && num == 1)
		std::cout << "row           :" << std::setw(15) << row << std::endl;
	else
	if (num == 1)
	{
		std::cout << "row           :" << std::setw(15) << row << ",\t";
		std::cout << "col           :" << std::setw(15) << col << std::endl;
	}
	else
	{
		std::cout << "row           :" << std::setw(15) << row << ",\t";
		std::cout << "col           :" << std::setw(15) << col << std::endl;
		std::cout << "num           :" << std::setw(15) << num << std::endl;
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
