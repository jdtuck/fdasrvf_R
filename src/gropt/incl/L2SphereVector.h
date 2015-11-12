#ifndef L2SPHEREVECTOR_H
#define L2SPHEREVECTOR_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

class L2SphereVector : public Element{
public:
	L2SphereVector(integer n);
	virtual L2SphereVector *ConstructEmpty(void) const;
};

#endif // end of STIEVECTOR_H
