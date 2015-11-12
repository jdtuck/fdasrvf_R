
#ifndef RBROYDENFAMILY_H
#define RBROYDENFAMILY_H

#include "SolversLS.h"
#include "def.h"

class RBroydenFamily : public SolversLS{
public:
	RBroydenFamily(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);
	virtual ~RBroydenFamily();
	double Phi(Variable *x2, Vector *y, Vector *s, LinearOPE *tildeH, double inpsy, double yHy, Vector *u);
	virtual void CheckParams();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool isconvex;
	double nu;
	double mu;
protected:
	bool isupdated;
	double betay, phic, inpsy, inpss;
	virtual void GetSearchDir();
	virtual void UpdateData();
	virtual void PrintInfo();

	Vector *s, *y, *u;
	LinearOPE *H, *tildeH;
};

#endif // end of RBROYDENFAMILY_H
