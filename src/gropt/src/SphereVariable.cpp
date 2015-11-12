
#include <SphereVariable.h>

SphereVariable::SphereVariable(integer n) :StieVariable(n)
{
};

SphereVariable *SphereVariable::ConstructEmpty(void) const
{
	return new SphereVariable(length);
};
