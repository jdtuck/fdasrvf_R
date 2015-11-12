
#include "OrthGroupVector.h"

OrthGroupVector::OrthGroupVector(integer n, integer m) :StieVector(n, m)
{
};

OrthGroupVector *OrthGroupVector::ConstructEmpty(void) const
{
	return new OrthGroupVector(size[0], size[1]);
};
