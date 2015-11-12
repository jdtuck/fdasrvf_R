
#ifndef RSD_H
#define RSD_H

#include "SolversLS.h"
#include "def.h"

class RSD : public SolversLS{
public:
	RSD(const Problem *prob, const Variable *initialx);
protected:
	virtual void GetSearchDir();
};

#endif // end of RSD_H
