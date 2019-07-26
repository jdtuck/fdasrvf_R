#include "ForDebug.h"

/*Define the namespace*/
namespace ROPTLIB{
	double ForDebug::NormF(const double *V, integer length)
	{
		return std::sqrt(ddot_(&length, const_cast<double *> (V), &GLOBAL::IONE, const_cast<double *> (V), &GLOBAL::IONE));
	};

	void ForDebug::Print(const char *name, const double *M, integer row, integer col, integer num)
	{
		Rprintf("=============%s============\n", name);
		if (col == 1 && num == 1)
		{
			for (integer i = 0; i < row; i++)
				Rprintf("%.10e\n", M[i]);
		}
		else
			if (num == 1)
			{
				for (integer j = 0; j < row; j++)
				{
					for (integer k = 0; k < col; k++)
					{
						Rprintf("%.10e\t", M[j + row * k]);
					}
					Rprintf("\n");
				}
			}
			else
			{
				for (integer i = 0; i < num; i++)
				{
					Rprintf("(:, :, %d)\n", i);
					for (integer j = 0; j < row; j++)
					{
						for (integer k = 0; k < col; k++)
						{
							Rprintf("%.10e\t", M[i * row * col + j + row * k]);
						}
						Rprintf("\n");
					}
				}
			}
	}
}; /*end of ROPTLIB namespace*/
