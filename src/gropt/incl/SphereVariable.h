#ifndef SPHEREVARIABLE_H
#define SPHEREVARIABLE_H

#include <StieVariable.h>

class SphereVariable : public StieVariable{
public:
	SphereVariable(integer n);
	virtual SphereVariable *ConstructEmpty(void) const;
};

#endif // end of SPHEREVARIABLE_H