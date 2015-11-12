
#ifndef SPLINE_H
#define SPLINE_H

#include <iostream>
#include <limits>
#include <cmath>
#include <def.h>

class Spline{
public:
	static int SplineUniformPeriodic(const double *Y, int n, double h, double *coefs); // periodic boundary condition: for closed curves
	static int SplinePeriodic(const double *X, const double *Y, int n, double *coefs);

	static int SplineUniformSlopes(const double *Y, int n, double h, double *coefs); // Hermite boundary conditions: for open curves (derivative-free from)
	static int SplineSlopes(const double *X, const double *Y, int n, double *coefs);

	static int SolveTridiagonalSystem(double *d, double *ud, double *ld, double *vec, double *s, int n);
	static int SolvePeriodicSystem(double *d, double *ud, double *ld, double *vec, double *s, int nn);

	static double ValSpline(const double *coefs, const double *breaks, int N, double t);
	static double ValSplineUniform(const double *coefs, int N, double h, double t);
	static void FirstDeri(const double *coefs, int N, double *dericoefs);
	static double ValFirstDeri(const double *dericoefs, int N, double h, double t);
	static void SecondDeri(const double *coefs, int N, double *dericoefs);
	static double ValSecondDeri(const double *dericoefs, int N, double h, double t);

};

#endif // SPLINE_H
