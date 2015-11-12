#ifndef OBLIQUEVARIABLE_H
#define OBLIQUEVARIABLE_H

#include <ProductElement.h>
#include <SphereVariable.h>

class ObliqueVariable : public ProductElement{
public:
	ObliqueVariable(integer n, integer num);
	virtual ~ObliqueVariable();
	virtual ObliqueVariable *ConstructEmpty(void) const;
	//virtual void RandInManifold();
	virtual void Print(const char *name = "", bool isonlymain = true) const;
};

#endif // end of OBLIQUEVARIABLE_H
