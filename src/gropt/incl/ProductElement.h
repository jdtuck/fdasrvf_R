
#ifndef PRODUCTELEMENT_H
#define PRODUCTELEMENT_H

#include "Element.h"
#include "Manifold.h"
#include <sstream>

class Element;

class ProductElement : public Element{
public:
	//void AddElement(Element *element, integer start, integer end);
	ProductElement(Element **elements, integer numofelements, integer *powsinterval, integer numoftypes);
	ProductElement(integer numberoftypes, ...);
	virtual ~ProductElement();

	virtual ProductElement *ConstructEmpty(void) const;
	virtual void CopyTo(Element *eta) const;
	virtual void RandUnform(double start = 0, double end = 1);
	virtual void RandGaussian(double mean = 0, double variance = 1);
	virtual void Print(const char *name = "", bool isonlymain = true) const;
	virtual void RandInManifold();

	virtual double *ObtainWriteEntireData();
	virtual double *ObtainWritePartialData();
	virtual void NewMemoryOnWrite();
	virtual void CopyOnWrite();
	virtual void CheckMemory(void) const;

	inline Element *GetElement(integer idx) const { return elements[idx]; };
	inline integer GetNumofElement(void) const { return numofelements; };

	//virtual void CopytoArray(double *array) const;
protected:
	virtual void ProductElementInitialization(Element **elements, integer numofelements, integer *powsinterval, integer numoftypes);
	virtual void ResetMemoryofElementsAndSpace();
	ProductElement();

	integer numoftypes;
	integer *powsinterval;

	Element **elements;
	integer numofelements;
};

#endif
