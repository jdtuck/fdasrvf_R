
#ifndef RTRNEWTON_H
#define RTRNEWTON_H

#include "SolversTR.h"
#include "def.h"

class RTRNewton : public SolversTR{
public:
	RTRNewton(const Problem *prob, const Variable *initialx);
protected:
	virtual void HessianEta(Vector *Eta, Vector *result);
};


#endif // end of RTRNEWTON_H