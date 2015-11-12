
#ifndef RWRBFGS_H
#define RWRBFGS_H

#include "SolversLS.h"
#include "def.h"

class RWRBFGS : public SolversLS{
public:
	RWRBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);
	virtual ~RWRBFGS();
	virtual void CheckParams();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool isconvex;
	double nu;
	double mu;
protected:
	bool isupdated;
	double inpsy, inpss;
	virtual void GetSearchDir();
	virtual void UpdateData();
	virtual void PrintInfo();

	Vector *s, *y;
	LinearOPE *H, *tildeH;
};

#endif // end of RWRBFGS_H
