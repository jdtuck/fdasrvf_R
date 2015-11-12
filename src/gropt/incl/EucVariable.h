#ifndef EUCVARIABLE_H
#define EUCVARIABLE_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

class EucVariable : public Element{
public:
	EucVariable(integer r, integer l = 1, integer n = 1);
	virtual EucVariable *ConstructEmpty(void) const;
	virtual void RandInManifold();
};

#endif // end of EUCVARIABLE_H
