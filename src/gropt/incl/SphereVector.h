#ifndef SPHEREVECTOR_H
#define SPHEREVECTOR_H

#include <StieVector.h>

class SphereVector :public StieVector{
public:
	SphereVector(integer n);
	virtual SphereVector *ConstructEmpty(void) const;
};

#endif // end of SPHEREVECTOR_H