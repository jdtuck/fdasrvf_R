
#ifndef DRIVERELASTICCURVESRO_H
#define DRIVERELASTICCURVESRO_H

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
#include <cmath>

#include <ElasticCurvesRO.h>

#include "def.h"
void DriverElasticCurvesRO(double *C1, double *C2, integer d, integer n, double w, bool rotated, bool isclosed,
	bool onlyDP, integer skipm, std::string solverstr, integer autoselectC, ProductElement *Xopt, bool &swap, double *fopts, double *comtime, integer &Nsout, integer &numinitialx);

double DynamicProgramming(const double *p_q1, const double *p_q2, integer d, integer N, double *gamma, bool isclosed);
void CenterC(double *C, integer d, integer n);
void NormalizedC(double *C, integer d, integer n);
double ComputeTotalAngle(const double *C, integer d, integer n);
void FindInitialBreaksAndNs(const double *C, integer d, integer n, integer minSkip, double thresholdsmall,
	integer rand_shift, integer *p_ms, integer &Lms, integer &Ns);
void CurveToQ(const double *C, integer d, integer n, double *q, bool isclosed);
void QToCurve(const double *Q, integer d, integer n, double *C, bool isclosed);
void ShiftC(const double *C, integer d, integer n, double *Cshift, integer m);
void FindBestRotation(const double *q1, const double *q2, integer d, integer n, double *O);
void GetCurveSmall(const double *C, double *Cs, integer d, integer n, integer ns, bool isclosed);
void GammaInverse(const double *DPgam, integer n, double *DPgamI);
void ReSampleGamma(const double *DPgams, integer ns, double *DPgam, integer n);
void Gradient(const double *DPgam, integer n, double h, double *grad);
void GradientPeriod(const double *DPgam, integer n, double h, double *grad);

#endif // end of DRIVERELASTICCURVESRO_H
