
#ifndef ELEMENT_H
#define ELEMENT_H

#include "randgen.h"
#include <cstdarg>
#include <map>
#include <string>
#include "SmartSpace.h"
#include "SharedSpace.h"
#include "def.h"

class SharedSpace;

typedef std::map<std::string, SharedSpace *> MAP;

class Element : public SmartSpace{
public:
	virtual ~Element();
	virtual Element *ConstructEmpty(void) const = 0;
	virtual void CopyTo(Element *eta) const;
	virtual void RandUnform(double start = 0, double end = 1);
	virtual void RandGaussian(double mean = 0, double variance = 1);
	virtual void Print(const char *name = "", bool isonlymain = true) const;
	virtual void RandInManifold();
	virtual double *ObtainWriteEntireData();
	virtual double *ObtainWritePartialData();

	virtual void AddToTempData(std::string name, SharedSpace * &Temp);
	virtual const SharedSpace *ObtainReadTempData(std::string name) const;
	virtual SharedSpace *ObtainWriteTempData(std::string name);
	virtual void RemoveFromTempData(std::string name);
	virtual void RemoveAllFromTempData();
	virtual bool TempDataExist(std::string name) const;
	//virtual void CopytoArray(double *array) const;
	void ObtainTempNames(std::string *names) const;

	inline integer GetSizeofTempData() const { return static_cast<integer> (TempData.size()); };
protected:
	MAP TempData;
};

#endif
