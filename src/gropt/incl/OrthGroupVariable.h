#ifndef ORTHGROUPVARIABLE_H
#define ORTHGROUPVARIABLE_H

#include <StieVariable.h>

class OrthGroupVariable : public StieVariable{
public:
	OrthGroupVariable(integer n);
	virtual OrthGroupVariable *ConstructEmpty(void) const;
};

#endif // end of ORTHGROUPVARIABLE_H