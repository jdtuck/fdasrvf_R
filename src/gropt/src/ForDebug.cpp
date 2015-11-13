
#include "ForDebug.h"

void ForDebug::Print(char *name, const double *M, integer row, integer col, integer num)
{
	Rcpp::Rcout << "=============" << name << "============" << std::endl;
	if (col == 1 && num == 1)
	{
		for (integer i = 0; i < row; i++)
			Rcpp::Rcout << M[i] << std::endl;
	}
	else
	if (num == 1)
	{
		for (integer j = 0; j < row; j++)
		{
			for (integer k = 0; k < col; k++)
			{
				Rcpp::Rcout << M[j + row * k] << "\t";
			}
			Rcpp::Rcout << std::endl;
		}
	}
	else
	{
		for (integer i = 0; i < num; i++)
		{
			Rcpp::Rcout << "(:, :, " << i << ")" << std::endl;
			for (integer j = 0; j < row; j++)
			{
				for (integer k = 0; k < col; k++)
				{
					Rcpp::Rcout << M[j + row * k] << "\t";
				}
				Rcpp::Rcout << std::endl;
			}
		}
	}
}
