#ifndef OBLIQUEVECTOR_H
#define OBLIQUEVECTOR_H

#include <ProductElement.h>
#include <SphereVector.h>

class ObliqueVector : public ProductElement{
public:
	ObliqueVector(integer n, integer num);
	virtual ~ObliqueVector();
	virtual ObliqueVector *ConstructEmpty(void) const;
	virtual void Print(const char *name = "", bool isonlymain = true) const;
};

#endif // end of OBLIQUEVECTOR_H
