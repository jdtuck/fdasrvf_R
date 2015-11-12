#ifndef ORTHGROUPVECTOR_H
#define ORTHGROUPVECTOR_H

#include <StieVector.h>

class OrthGroupVector :public StieVector{
public:
	OrthGroupVector(integer n, integer m = 1);
	virtual OrthGroupVector *ConstructEmpty(void) const;
};

#endif // end of ORTHGROUPVECTOR_H
