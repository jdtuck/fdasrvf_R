
#ifndef LRTRSR1_H
#define LRTRSR1_H

#include "ForDebug.h"
#include "SolversTR.h"
#include "def.h"

class LRTRSR1 : public SolversTR{
public:
	LRTRSR1(const Problem *prob, const Variable *initialx);
	virtual ~LRTRSR1();
	virtual void CheckParams();
	virtual void Run();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	bool isconvex;
	integer LengthSY;
protected:
	virtual void HessianEta(Vector *Eta, Vector *result);
	virtual void UpdateData();
	virtual void Acceptence();
	virtual void PrintInfo();

	bool isupdated;
	bool ischangedSandY;
	double inpsy, inpss, inpyy;
	Vector *s, *y;
	Vector **S, **Y, **YMGS;
	double *SS, *SY, *PMGQ;
	integer *P;
	double gamma;
	integer Currentlength;
	integer beginidx;
private:
	void NewVectors(Vector ** &Vs, integer l);
	void DeleteVectors(Vector ** &Vs, integer l);
};

#endif // end of LRTRSR1_H
