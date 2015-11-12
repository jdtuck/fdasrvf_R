#ifndef SHAREDSPACE_H
#define SHAREDSPACE_H

#include "SmartSpace.h"
#include <Element.h>
#include "def.h"

class Element;

class SharedSpace : public SmartSpace{
public:
	SharedSpace(integer numberofdimensions, ...);
	SharedSpace(Element *inelement);
	~SharedSpace();
	virtual SharedSpace *ConstructEmpty(void) const;
	virtual void CopyTo(SharedSpace *eta) const;
	virtual void Print(const char *name = "") const;
	Element *GetSharedElement(void) const;
private:
	Element *SharedElement;
};

#endif // end of SHAREDSPACE_H
