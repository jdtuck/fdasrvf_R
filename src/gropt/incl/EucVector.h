#ifndef EUCVECTOR_H
#define EUCVECTOR_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

class EucVector : public Element{
public:
	EucVector(integer r, integer c = 1, integer n = 1);
	virtual EucVector *ConstructEmpty(void) const;
};


#endif // end of EUCVECTOR_H
