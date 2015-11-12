
#include "ForDebug.h"

void ForDebug::Print(char *name, const double *M, integer row, integer col, integer num)
{
	std::cout << "=============" << name << "============" << std::endl;
	if (col == 1 && num == 1)
	{
		for (integer i = 0; i < row; i++)
			std::cout << M[i] << std::endl;
	}
	else
	if (num == 1)
	{
		for (integer j = 0; j < row; j++)
		{
			for (integer k = 0; k < col; k++)
			{
				std::cout << M[j + row * k] << "\t";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		for (integer i = 0; i < num; i++)
		{
			std::cout << "(:, :, " << i << ")" << std::endl;
			for (integer j = 0; j < row; j++)
			{
				for (integer k = 0; k < col; k++)
				{
					std::cout << M[j + row * k] << "\t";
				}
				std::cout << std::endl;
			}
		}
	}
}
