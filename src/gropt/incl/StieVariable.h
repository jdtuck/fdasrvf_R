#ifndef STIEVARIABLE_H
#define STIEVARIABLE_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

class StieVariable : public Element{
public:
	StieVariable(integer r, integer l = 1, integer n = 1);
	virtual StieVariable *ConstructEmpty(void) const;
	virtual void RandInManifold();
};

#endif // end of STIEVARIABLE_H
