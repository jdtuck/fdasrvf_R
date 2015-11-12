
#ifndef OBLIQUE_H
#define OBLIQUE_H

#include <ProductManifold.h>
#include <Sphere.h>
#include "ObliqueVariable.h"
#include "ObliqueVector.h"

// It is the products of unit spheres
class Oblique : public ProductManifold{
public:
	Oblique(integer n, integer num);
	~Oblique();

	// Default one: choose qf, parallelization and intrinsic approach and use householder reflections
	virtual void ChooseObliqueParamsSet1();

	// choose exponential map, parallel translation and extrinsic approach and no householder reflections
	// Even though the Householder reflections are not used, the locking condition is satisfied.
	void ChooseObliqueParamsSet2();

	// choose qf, parallel translation and extrinsic approach and no householder reflections
	// The locking conidition is not satisfied
	void ChooseObliqueParamsSet3();

	// choose qf, parallel translation and extrinsic approach and no householder reflections
	// Beta \neq 1 is used and the locking conidition is satisfied
	void ChooseObliqueParamsSet4();

	virtual void SetParams(PARAMSMAP params);
};

#endif // end of OBLIQUE_H
