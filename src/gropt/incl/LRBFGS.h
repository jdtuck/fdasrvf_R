
#ifndef LRBFGS_H
#define LRBFGS_H

#include <cstring>
#include "SolversLS.h"
#include "def.h"

class LRBFGS : public SolversLS{
public:
	LRBFGS(const Problem *prob, const Variable *initialx);
	virtual ~LRBFGS();
	virtual void CheckParams();
	virtual void Run();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool isconvex;
	double nu;
	double mu;
	integer LengthSY;
protected:
	bool isupdated;
	double betay, inpsy, inpss, inpyy;
	virtual void GetSearchDir();
	virtual void UpdateData();
	virtual void PrintInfo();

	Vector *s, *y;
	Vector **S, **Y;
	double *RHO;
	double rho, gamma;
	integer Currentlength;
	integer beginidx;
private:
	void NewVectors(Vector ** &Vs, integer l);
	void DeleteVectors(Vector ** &Vs, integer l);
};

#endif // end of RBROYDENFAMILY_H
