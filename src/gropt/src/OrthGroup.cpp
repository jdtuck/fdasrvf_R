
#include "OrthGroup.h"

OrthGroup::OrthGroup(integer inn) :Stiefel(inn, inn)
{
	name.assign("OrthGroup");
	delete EMPTYEXTR;
	delete EMPTYINTR;
	EMPTYEXTR = new OrthGroupVector(n, n);
//	std::cout << "n:" << n << std::endl;//---
//	std::cout << "IntrinsicDim:" << IntrinsicDim << std::endl;//---
	EMPTYINTR = new OrthGroupVector(IntrinsicDim);
};

OrthGroup::~OrthGroup(void)
{
};
