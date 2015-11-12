
#include "LowRankVector.h"

LowRankVector::LowRankVector(integer Ur, integer Uc, integer Drc, integer Vr, integer Vc)
{
	StieVector dU(Ur, Uc);
	EucVector dD(Drc, Drc);
	StieVector dV(Vr, Vc);

	Element **Elems = new Element *[3];
	Elems[0] = &dU;
	Elems[1] = &dD;
	Elems[2] = &dV;
	integer *powsintev = new integer[4];
	powsintev[0] = 0;
	powsintev[1] = 1;
	powsintev[2] = 2;
	powsintev[3] = 3;

	ProductElementInitialization(Elems, 3, powsintev, 3);

	delete[] powsintev;
	delete[] Elems;
};

LowRankVector::~LowRankVector(void)
{
};

LowRankVector *LowRankVector::ConstructEmpty(void) const
{
	integer Ur = elements[0]->Getsize()[0];
	integer Uc = elements[0]->Getsize()[1];
	integer Drc = elements[1]->Getsize()[0];
	integer Vr = elements[2]->Getsize()[0];
	integer Vc = elements[2]->Getsize()[1];
	return new LowRankVector(Ur, Uc, Drc, Vr, Vc);
};
