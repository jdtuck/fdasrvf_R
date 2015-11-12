#ifndef L2SPHEREVARIABLE_H
#define L2SPHEREVARIABLE_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

class L2SphereVariable : public Element{
public:
	L2SphereVariable(integer n);
	virtual L2SphereVariable *ConstructEmpty(void) const;
	virtual void RandInManifold();
};

#endif // end of L2SPHEREVARIABLE_H
