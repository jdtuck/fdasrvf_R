#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Manifold.h"
#include "Problem.h"
#include "def.h"

enum StopCrit{ FUN_REL, GRAD_F, GRAD_F_0, STOPCRITLENGTH };
enum DEBUGINFO{ NOOUTPUT, FINALRESULT, ITERRESULT, DETAILED, DEBUGLENGTH };

class Solvers{
public:
	virtual void Run();
	virtual void CheckParams();
	virtual void OutPutResults(Variable *inx1, double &inf1, double &inngf0, double &inngf, integer &initer,
		integer &innf, integer &inng, integer &innR, integer &innV, integer &innVp, double &inComTime,
		double *intimeSeries, double *infunSeries, double *ingradSeries, integer &inlengthSeries);

	inline const Variable *GetXopt(void) const { return x1; };
	inline double Getfinalfun(void) const { return f2; };
	inline double Getnormgf(void) const { return ngf; };
	inline double Getnormgfgf0(void) const { return ngf / ngf0; };
	inline double GetComTime(void) const { return ComTime; };
	inline integer Getnf(void) const { return nf; };
	inline integer Getng(void) const { return ng; };
	inline integer GetnR(void) const { return nR; };
	inline integer GetnV(void) const { return nV; };
	inline integer GetnVp(void) const { return nVp; };
	inline integer GetnH(void) const { return nH; };
	inline integer GetIter(void) const { return iter; };
	inline integer GetlengthSeries(void) const { return lengthSeries; };
	inline double *GettimeSeries(void) const { return timeSeries; };
	inline double *GetfunSeries(void) const { return funSeries; };
	inline double *GetgradSeries(void) const { return gradSeries; };

	virtual void SetParams(PARAMSMAP params);
	// parameters
	// Solvers
	StopCrit Stop_Criterion;
	double Tolerance;
	integer Max_Iteration;
	integer Min_Iteration;
	integer OutputGap;
	DEBUGINFO DEBUG;

	virtual ~Solvers(void) = 0;
protected:
	virtual void Initialization(const Problem *prob, const Variable *initialx, const Vector *EMPTYETA);
	virtual void UpdateData(void) = 0;
	virtual bool IsStopped(void);
	virtual void PrintGenInfo(void);
	virtual void PrintInfo(void);

	// algorithm-related variables:
	Variable *x1, *x2;
	Vector *gf1, *gf2;
	double f1, f2;
	double ngf0, ngf;

	// Input parameters and functions
	const Manifold *Mani;
	const Problem *Prob;

	// For debug information
	integer iter;
	unsigned long starttime;
	double ComTime;
	integer nf, ng, nR, nV, nVp, nH;
	double *timeSeries, *funSeries, *gradSeries;
	integer lengthSeries;
	std::string SolverName;

	integer numfref;
};

#endif // end of SOLVER_H
