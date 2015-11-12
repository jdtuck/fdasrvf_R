#ifndef SOLVERLS_H
#define SOLVERLS_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Manifold.h"
#include "Problem.h"
#include "Solvers.h"
#include "def.h"

enum LSAlgo{ ARMIJO, WOLFE, STRONGWOLFE, EXACT, LSALGOLENGTH };
enum LSstatusSet{NOCURVATURE, MINSTEPSIZE, MAXSTEPSIZE, NONEXACT, LSERROR, SUCCESS, LSSTATUSSETLENGTH};

class SolversLS : public Solvers{
public:
	virtual void Run();
	virtual void CheckParams();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	LSAlgo LineSearch_LS;
	double LS_alpha;
	double LS_beta;
	double Minstepsize;
	double Maxstepsize;
	double LS_ratio1;
	double LS_ratio2;
	double Initstepsize;
protected:
	virtual ~SolversLS();
	virtual void GetSearchDir(void) = 0; // required to be overload in derived class

	virtual void Initialization(const Problem *prob, const Variable *initialx, const Vector *EMPTYETA);
	virtual void LinesearchArmijo();
	virtual void LinesearchWolfe();
	virtual void LinesearchStrongWolfe();
	virtual void LinesearchExact();
	virtual double h();
	virtual double dh();
	virtual void UpdateData();

	void (SolversLS::*Linesearch)();

	// algorithm-related variables:
	Vector *eta1, *eta2;
	Vector *zeta; // temp storage

	// parameters
	double initiallength;
	double stepsize;
	double initialslope, newslope;
	double fpre;
	LSstatusSet LSstatus;
	std::string *LSstatusSetnames;
private:
	virtual void PrintGenInfo();
	void Zoom(double x1, double fx1, double slopex1, double x2, double fx2);
};

#endif // end of SOLVERLS_H
