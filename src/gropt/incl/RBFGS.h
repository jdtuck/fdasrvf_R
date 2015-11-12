
#ifndef RBFGS_H
#define RBFGS_H

#include "SolversLS.h"
#include "def.h"

class RBFGS : public SolversLS{
public:
	RBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);
	virtual ~RBFGS();
	virtual void CheckParams();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool isconvex;
	double nu;
	double mu;
protected:
	bool isupdated;
	double betay, inpsy, inpss;
	virtual void GetSearchDir();
	virtual void UpdateData();
	virtual void PrintInfo();

	Vector *s, *y;
	LinearOPE *H, *tildeH;
};

#endif // end of RBFGS_H
