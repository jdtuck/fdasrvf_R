#ifndef STIEVECTOR_H
#define STIEVECTOR_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

class StieVector : public Element{
public:
	StieVector(integer r, integer c = 1, integer n = 1);
	virtual StieVector *ConstructEmpty(void) const;
};


#endif // end of STIEVECTOR_H
