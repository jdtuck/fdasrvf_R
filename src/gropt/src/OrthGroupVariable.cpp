
#include "OrthGroupVariable.h"

OrthGroupVariable::OrthGroupVariable(integer n) :StieVariable(n, n)
{
};

OrthGroupVariable *OrthGroupVariable::ConstructEmpty(void) const
{
	return new OrthGroupVariable(size[0]);
};
