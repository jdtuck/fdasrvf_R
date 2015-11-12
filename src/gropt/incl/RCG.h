
#ifndef RCG_H
#define RCG_H

#include "SolversLS.h"
#include "def.h"
#undef max

enum RCGmethods{ FLETCHER_REEVES, POLAK_RIBIERE_MOD, HESTENES_STIEFEL, FR_PR, DAI_YUAN, HAGER_ZHANG, RCGMETHODSLENGTH };

class RCG : public SolversLS{
public:
	RCG(const Problem *prob, const Variable *initialx);
	virtual ~RCG();
	virtual void CheckParams();

	virtual void SetParams(PARAMSMAP params);
	// parameters
	integer ManDim;
	RCGmethods RCGmethod;
protected:
	virtual void GetSearchDir();
	virtual void UpdateData();
	virtual void PrintInfo();
	std::string *RCGmethodSetnames;
	double sigma;
};

#endif // end of RCG_H
