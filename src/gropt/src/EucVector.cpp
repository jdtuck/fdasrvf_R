
#include "EucVector.h"

EucVector::EucVector(integer r, integer c, integer n)
{
	Element::Initialization(3, r, c, n);
};

EucVector *EucVector::ConstructEmpty(void) const
{
	return new EucVector(size[0], size[1], size[2]);
};
