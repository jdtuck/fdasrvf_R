#ifndef DRIVERMEXPROB_H
#define DRIVERMEXPROB_H

#include "ForDebug.h"
#include <iostream>
#include "randgen.h"
#include "Manifold.h"
#include "Problem.h"
#include "SolversLS.h"
#include <ctime>

#include "EucVariable.h"
#include "EucVector.h"
#include "EucFrechetMean.h"
#include "EucQuadratic.h"

#include "StieBrockett.h"
#include "StieVector.h"
#include "StieVariable.h"
#include "Stiefel.h"

#include "RSD.h"
#include "RNewton.h"
#include "RCG.h"
#include "RBroydenFamily.h"
#include "RWRBFGS.h"
#include "RBFGS.h"
#include "LRBFGS.h"

#include "SolversTR.h"
#include "RTRSD.h"
#include "RTRNewton.h"
#include "RTRSR1.h"
#include "LRTRSR1.h"

#include "def.h"

#include <Sphere.h>
#include <SphereVariable.h>
#include <SphereVector.h>

#include <Oblique.h>
#include <ObliqueVariable.h>
#include <ObliqueVector.h>

#include <LowRank.h>
#include <LowRankVariable.h>
#include <LowRankVector.h>

#include <ProductElement.h>
#include <ProductManifold.h>

#ifdef MATLAB_MEX_FILE

#include <mexProblem.h>

void DriverMexProb(int &nlhs, mxArray ** &plhs, int &nrhs, const mxArray ** &prhs);
void ParseSolverParamsAndOptimizing(const mxArray *SolverParams, Problem *Prob, Variable *initialX, mxArray ** &plhs);
bool ParseManiParams(const mxArray *ManiParams, Manifold **&manifolds, Element **&elements,
	integer *&powsinterval, integer &numoftype, integer &numoftotal);
Manifold *GetAManifold(const char *name, integer n, integer m, integer p = 1);
Element *GetAnElement(const char *name, integer n, integer m, integer p = 1);

void ParseSolverParamsAndOptimizing(const mxArray *SolverParams, Problem *Prob, Variable *initialX, mxArray **&plhs)
{
	integer nfields = mxGetNumberOfFields(SolverParams);
	const char *name;
	mxArray *tmp;
	double value;
	PARAMSMAP params;
	std::string key;
	for (integer i = 0; i < nfields; i++)
	{
		name = mxGetFieldNameByNumber(SolverParams, i);
		tmp = mxGetFieldByNumber(SolverParams, 0, i);
		value = mxGetScalar(tmp);
		key.assign(name);
		params.insert(std::pair<std::string, double>(key, value));
	}

	tmp = mexProblem::GetFieldbyName(SolverParams, 0, "method");
	if (tmp == nullptr)
	{
		mexErrMsgTxt("A method must be specified.");
	}
	char methodname[30] = "";
	mxGetString(tmp, methodname, 30);
	std::string stdmethodname = methodname;
	Solvers *solver;
	if (stdmethodname == "RSD")
	{
		solver = new RSD(Prob, initialX);
	}
	else
	if (stdmethodname == "RNewton")
	{
		solver = new RNewton(Prob, initialX);
	}
	else
	if (stdmethodname == "RCG")
	{
		solver = new RCG(Prob, initialX);
	}
	else
	if (stdmethodname == "RBroydenFamily")
	{
		solver = new RBroydenFamily(Prob, initialX);
	}
	else
	if (stdmethodname == "RWRBFGS")
	{
		solver = new RWRBFGS(Prob, initialX);
	}
	else
	if (stdmethodname == "RBFGS")
	{
		solver = new RBFGS(Prob, initialX);
	}
	else
	if (stdmethodname == "LRBFGS")
	{
		solver = new LRBFGS(Prob, initialX);
	}
	else
	if (stdmethodname == "RTRSD")
	{
		solver = new RTRSD(Prob, initialX);
	}
	else
	if (stdmethodname == "RTRNewton")
	{
		solver = new RTRNewton(Prob, initialX);
	}
	else
	if (stdmethodname == "RTRSR1")
	{
		solver = new RTRSR1(Prob, initialX);
	}
	else
	if (stdmethodname == "LRTRSR1")
	{
		solver = new LRTRSR1(Prob, initialX);
	}
	else
	{
		std::cout << "Warning: Unrecognized solver: " << stdmethodname << ". Use LRBFGS instead!" << std::endl;
		solver = new LRBFGS(Prob, initialX);
	}
	solver->SetParams(params);

	tmp = mexProblem::GetFieldbyName(SolverParams, 0, "IsCheckParams");
	if (tmp != nullptr)
	{
		if (fabs(mxGetScalar(tmp)) > std::numeric_limits<double>::epsilon()) // if the value is nonzero
		{
			solver->CheckParams();
		}
	}
	solver->Run();
	mexProblem::ObtainMxArrayFromElement(plhs[0], solver->GetXopt());
	plhs[1] = mxCreateDoubleScalar(static_cast<double> (solver->Getfinalfun()));
	plhs[2] = mxCreateDoubleScalar(static_cast<double> (solver->Getnormgf()));
	plhs[3] = mxCreateDoubleScalar(static_cast<double> (solver->Getnormgfgf0()));
	plhs[4] = mxCreateDoubleScalar(static_cast<double> (solver->GetIter()));
	plhs[5] = mxCreateDoubleScalar(static_cast<double> (solver->Getnf()));
	plhs[6] = mxCreateDoubleScalar(static_cast<double> (solver->Getng()));
	plhs[7] = mxCreateDoubleScalar(static_cast<double> (solver->GetnR()));
	plhs[8] = mxCreateDoubleScalar(static_cast<double> (solver->GetnV()));
	plhs[9] = mxCreateDoubleScalar(static_cast<double> (solver->GetnVp()));
	plhs[10] = mxCreateDoubleScalar(static_cast<double> (solver->GetnH()));
	plhs[11] = mxCreateDoubleScalar(static_cast<double> (solver->GetComTime()));
	integer lengthSeries = solver->GetlengthSeries();
	plhs[12] = mxCreateDoubleMatrix(lengthSeries, 1, mxREAL);
	plhs[13] = mxCreateDoubleMatrix(lengthSeries, 1, mxREAL);
	plhs[14] = mxCreateDoubleMatrix(lengthSeries, 1, mxREAL);
	double *plhsfun = mxGetPr(plhs[12]), *plhsgrad = mxGetPr(plhs[13]), *plhstime = mxGetPr(plhs[14]);
	for (integer i = 0; i < lengthSeries; i++)
	{
		plhsfun[i] = solver->GetfunSeries()[i];
		plhsgrad[i] = solver->GetgradSeries()[i];
		plhstime[i] = solver->GettimeSeries()[i];
	}

	tmp = mexProblem::GetFieldbyName(SolverParams, 0, "IsCheckGradHess");
	if (tmp != nullptr)
	{
		if (fabs(mxGetScalar(tmp)) > std::numeric_limits<double>::epsilon()) // if the value is nonzero
		{
			// Check gradient and Hessian
			const Variable *xopt = solver->GetXopt();
			Prob->CheckGradHessian(initialX);
			Prob->CheckGradHessian(xopt);
		}
	}

	delete solver;
};

bool ParseManiParams(const mxArray *ManiParams, Manifold **&manifolds, Element **&elements,
	integer *&powsinterval, integer &numoftype, integer &numoftotal)
{
	// Parse ManiParams
	integer nfields4 = mxGetNumberOfFields(ManiParams);
	numoftype = mxGetNumberOfElements(ManiParams);
	powsinterval = new integer[numoftype + 1];
	char name[30] = "";
	manifolds = new Manifold *[numoftype];
	integer n, p, m, Params;
	powsinterval[0] = 0;

	for (integer i = 0; i < numoftype; i++)
		powsinterval[i + 1] = powsinterval[i] + mxGetScalar(mexProblem::GetFieldbyName(ManiParams, i, "numofmani"));

	numoftotal = powsinterval[numoftype];
	elements = new Element *[numoftotal];
	PARAMSMAP params;
	for (integer i = 0; i < numoftype; i++)
	{
		if (mxGetString(mexProblem::GetFieldbyName(ManiParams, i, "name"), name, 30))
			mexErrMsgTxt("error in getting manifold name!");
		n = mxGetScalar(mexProblem::GetFieldbyName(ManiParams, i, "n"));
		m = mxGetScalar(mexProblem::GetFieldbyName(ManiParams, i, "m"));
		p = mxGetScalar(mexProblem::GetFieldbyName(ManiParams, i, "p"));
		Params = mxGetScalar(mexProblem::GetFieldbyName(ManiParams, i, "ParamSet"));
        params[static_cast<std::string> ("ParamSet")] = Params;
//		params.insert(std::pair<std::string, double>("ParamSet", Params));

		manifolds[i] = GetAManifold(name, n, m, p);
		manifolds[i]->SetParams(params);

		for (integer j = powsinterval[i]; j < powsinterval[i + 1]; j++)
		{
			elements[j] = GetAnElement(name, n, m, p);
		}
		if (manifolds[i] == nullptr || elements[i] == nullptr)
		{
			return false;
		}
	}

	return true;
};

Manifold *GetAManifold(const char *name, integer n, integer m, integer p)
{
	if (strcmp(name, "Euclidean") == 0)
	{
		return new Euclidean(n, m);
	}
	else
	if (strcmp(name, "Sphere") == 0)
	{
		return new Sphere(n);
	}
	else
	if (strcmp(name, "Stiefel") == 0)
	{
		return new Stiefel(n, p);
	}
	else
	if (strcmp(name, "Oblique") == 0)
	{
		return new Oblique(n, m);
	}
	else
	if (strcmp(name, "LowRank") == 0)
	{
		return new LowRank(n, m, p);
	}
	else
	{
		std::cout << "Manifold: " << name << " does not implemented in this library!" << std::endl;
		return nullptr;
	}
};

Element *GetAnElement(const char *name, integer n, integer m, integer p)
{
	if (strcmp(name, "Euclidean") == 0)
	{
		return new EucVariable(n, m);
	}
	else
	if (strcmp(name, "Sphere") == 0)
	{
		return new SphereVariable(n);
	}
	else
	if (strcmp(name, "Stiefel") == 0)
	{
		return new StieVariable(n, p);
	}
	else
	if (strcmp(name, "Oblique") == 0)
	{
		return new ObliqueVariable(n, m);
	}
	else
	if (strcmp(name, "LowRank") == 0)
	{
		return new LowRankVariable(n, m, p);
	}
	else
	{
		std::cout << "Element: " << name << " does not implemented in this library!" << std::endl;
		return nullptr;
	}
};

#endif // end of MATLAB_MEX_FILE

#endif // end of DRIVERMEXPROB_H
