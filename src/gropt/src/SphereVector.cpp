
#include <SphereVector.h>

SphereVector::SphereVector(integer n) :StieVector(n)
{
};

SphereVector *SphereVector::ConstructEmpty(void) const
{
	return new SphereVector(length);
};
