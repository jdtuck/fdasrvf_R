#ifndef TESTTESTSPARSEPCA_H
#define TESTTESTSPARSEPCA_H


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

#include <ObliqueTestSparsePCA.h>
#include <Oblique.h>
#include <ObliqueVariable.h>
#include <ObliqueVector.h>

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

integer GetNumBetweenC1andC2(const Element *x, double c1, double c2);
void testTestSparsePCA(double *B, double *Dsq, integer p, integer n, integer r, double epsilon, double mu, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTTESTSPARSEPCA_H