#ifndef FORDEBUG_H
#define FORDEBUG_H

#include <iostream>
#include "def.h"

class ForDebug{
public:
	static void Print(char *name, const double *M, integer row, integer col = 1, integer num = 1);
};

#endif // end of FORDEBUG_H
