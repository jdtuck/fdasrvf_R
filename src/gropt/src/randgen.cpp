
#ifndef RANDGEN_CPP
#define RANDGEN_CPP

#include "randgen.h"
#include "def.h"

void genrandseed(unsigned int s) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(s));
}

double genrandreal(void)
{
    Rcpp::NumericVector X(1);
    X = Rcpp::runif(1);
    return static_cast<double> (X[0]) / RAND_MAX;
}

double genrandnormal(void)
{
	static double rand1, rand2;
	double tmp = genrandreal();
	while (tmp == 1.0)
		tmp = genrandreal();
    rand1 = -2 * log(1.0 - tmp);
    rand2 = (1.0 - genrandreal()) * 6.2831853071795864769252866;
	return sqrt(rand1) * cos(rand2);
}

#endif
