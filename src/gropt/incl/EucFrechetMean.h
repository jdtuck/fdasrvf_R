
#ifndef EUCFRECHETMEAN_H
#define EUCFRECHETMEAN_H

#include "Euclidean.h"
#include "EucVariable.h"
#include "EucVector.h"
#include "Problem.h"
#include "def.h"

// problem: \min_{x \in R^dim} \sum_{i = 1}^num W(i) \|x - D(:, i)\|_2^2
class EucFrechetMean : public Problem{
public:
	EucFrechetMean(double *W, double *D, integer num, integer dim);
	virtual ~EucFrechetMean();
	virtual double f(Variable *x) const;
	virtual void Grad(Variable *x, Vector *gf) const;
	virtual void HessianEta(Variable *x, Vector *etax, Vector *xix) const;

	double *Weights; // length num
	double *Data;    // length num * dim
	integer Num;
	integer Dim;
};

#endif
