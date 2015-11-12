
#ifndef SMARTSPACE_H
#define SMARTSPACE_H

#include "randgen.h"
#include <cstdarg>
#include <map>
#include "def.h"

#ifdef CHECKMEMORYDELETED
extern std::map<integer *, integer> *CheckMemoryDeleted;
#endif

class SmartSpace{
public:
	virtual SmartSpace *ConstructEmpty(void) const = 0;
	virtual void Initialization(integer numberofdimensions, ...);
	virtual void CopyTo(SmartSpace *eta) const;
	virtual void RandUnform(double start = 0, double end = 1);
	virtual void RandGaussian(double mean = 0, double variance = 1);
	virtual const double *ObtainReadData(void) const;
	virtual double *ObtainWriteEntireData(void);
	virtual double *ObtainWritePartialData(void);
	virtual void NewMemoryOnWrite(void);
	virtual void CopyOnWrite(void);

	inline const integer *GetSharedTimes(void) const { return sharedtimes; };
	inline const double *GetSpace(void) const { return Space; };

	virtual ~SmartSpace(void) = 0;
	virtual void Print(const char *name = "") const;

	inline const integer *Getsize(void) const { return size; };
	inline integer Getls(void) const { return ls; };
	inline integer Getlength(void) const { return length; };

	virtual void SetByParams(integer *size, integer ls, integer length, integer *sharedtimes, double *Space);
	virtual void DeleteBySettingNull(void);
protected:
	integer *size;
	integer ls; // length of size
	integer length;

	integer *sharedtimes;
	double *Space;

	void NewMemory(void);
};

#endif
