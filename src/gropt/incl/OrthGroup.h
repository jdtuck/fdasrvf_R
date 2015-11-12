#ifndef ORTHGROUP_H
#define ORTHGROUP_H

#include "OrthGroupVariable.h"
#include "OrthGroupVector.h"
#include <Stiefel.h>

class OrthGroup : public Stiefel{
public:
	OrthGroup(integer inn);
	virtual ~OrthGroup();
};

#endif // end of ORTHGROUP_H
