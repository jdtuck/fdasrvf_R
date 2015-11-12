#ifndef LOWRANKVECTOR_H
#define LOWRANKVECTOR_H

#include <ProductElement.h>
#include <StieVector.h>
#include <EucVector.h>

class LowRankVector : public ProductElement{
public:
	LowRankVector(integer Ur, integer Uc, integer Drc, integer Vr, integer Vc);
	virtual ~LowRankVector(void);
	virtual LowRankVector *ConstructEmpty(void) const;
};

#endif // end of OBLIQUEVECTOR_H
