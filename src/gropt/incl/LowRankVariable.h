#ifndef LOWRANKVARIABLE_H
#define LOWRANKVARIABLE_H

#include <ProductElement.h>
#include <StieVariable.h>
#include <EucVariable.h>

class LowRankVariable : public ProductElement{
public:
	LowRankVariable(integer m, integer n, integer r);
	virtual ~LowRankVariable(void);
	virtual LowRankVariable *ConstructEmpty(void) const;
};

#endif // end of OBLIQUEVARIABLE_H
