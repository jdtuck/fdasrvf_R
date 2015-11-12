
#ifndef RTRSD_H
#define RTRSD_H

#include "SolversTR.h"
#include "def.h"

class RTRSD : public SolversTR{
public:
	RTRSD(const Problem *prob, const Variable *initialx);
protected:
	virtual void HessianEta(Vector *Eta, Vector *result);
};


#endif // end of RTRSD_H
