
#ifndef LINEAROPE_H
#define LINEAROPE_H

#include "Element.h"
#include "def.h"

class LinearOPE : public SmartSpace{
public:
	LinearOPE(integer s);
	virtual LinearOPE *ConstructEmpty(void) const;
	virtual void ScaledIdOPE(double scaler = 1);
};

#endif // end of LINEAROPE_H
