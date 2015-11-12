
#ifndef SOLVERSTR_H
#define SOLVERSTR_H

#include "Solvers.h"
#include "def.h"

enum tCGstatusSet{TR_NEGCURVTURE, TR_EXCREGION, TR_LCON, TR_SCON, TR_MAXITER, TCGSTATUSSETLENGTH};

class SolversTR : public Solvers{
public:
	virtual void CheckParams();
	virtual void Run();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	double Acceptence_Rho;
	double Shrinked_tau;
	double Magnified_tau;
	double minimum_Delta;
	double maximum_Delta;
	bool useRand;
	integer Max_Inner_Iter;
	integer Min_Inner_Iter;
	double theta;
	double kappa;
	double initial_Delta;
protected:
	virtual ~SolversTR();
	virtual void HessianEta(Vector *Eta, Vector *result) = 0; // required to be overloaded in derived class

	virtual void InitialVector();
	virtual void tCG_TR();
	virtual void PreConditioner(Variable *x, Vector *eta, Vector *result);
	virtual void Initialization(const Problem *prob, const Variable *initialx, const Vector *EMPTYETA);
	virtual void UpdateData();
	virtual void Acceptence();

	// algorithm-related variables:
	Vector *eta1, *eta2;
	Vector *zeta; // temp storage
	Vector *r, *z, *delta, *Hd;
	double rho;
	double Delta;
	integer innerIter;
	tCGstatusSet tCGstatus;
	std::string *tCGstatusSetnames;
private:
	virtual void PrintGenInfo();
};

#endif
