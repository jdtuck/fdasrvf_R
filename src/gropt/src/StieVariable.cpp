
#include "StieVariable.h"

StieVariable::StieVariable(integer r, integer c, integer n)
{
	Element::Initialization(3, r, c, n);
};

StieVariable *StieVariable::ConstructEmpty(void) const
{
	return new StieVariable(size[0], size[1], size[2]);
};

void StieVariable::RandInManifold(void)
{
	this->RandGaussian();
	double *xU = this->ObtainWriteEntireData();
	integer n = size[0], p = size[1];
	integer *jpvt = new integer[p];
	integer lwork = 2 * p + (1 + p) * INITIALBLOCKSIZE, info;
	double *tau = new double[p + lwork];
	double *work = tau + p;
	for (integer i = 0; i < p; i++)
		jpvt[i] = 0;
	dgeqp3_(&n, &p, xU, &n, jpvt, tau, work, &lwork, &info);
	if (info < 0)
		std::cout << "Error in qr decomposition!" << std::endl;
	dorgqr_(&n, &p, &p, xU, &n, tau, work, &lwork, &info);
	if (info < 0)
		std::cout << "Error in forming Q matrix!" << std::endl;
	delete[] jpvt;
	delete[] tau;
};
