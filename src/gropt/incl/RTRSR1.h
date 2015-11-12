
#ifndef RTRSR1_H
#define RTRSR1_H

#include "SolversTR.h"
#include "def.h"

class RTRSR1 : public SolversTR{
public:
	RTRSR1(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr);
	virtual ~RTRSR1();
	virtual void CheckParams();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool isconvex;
protected:
	virtual void HessianEta(Vector *Eta, Vector *result);
	virtual void UpdateData();
	virtual void Acceptence();
	virtual void PrintInfo();

	double inpss;
	bool isupdated;
	Vector *s, *y;
	LinearOPE *B, *tildeB;
};


#endif // end of RTRSR1_H
