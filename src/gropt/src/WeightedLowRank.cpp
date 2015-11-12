
#include "WeightedLowRank.h"

WeightedLowRank::WeightedLowRank(double *inA, double *inW, integer inm, integer inn, integer inr)
{
	A = inA;
	W = inW;
	m = inm;
	n = inn;
	r = inr;
};

WeightedLowRank::~WeightedLowRank(void)
{
};

double WeightedLowRank::f(Variable *x) const
{
    ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
    const double *Uptr = ProdxxM->GetElement(0)->ObtainReadData();
   	const double *Dptr = ProdxxM->GetElement(1)->ObtainReadData();
    const double *Vptr = ProdxxM->GetElement(2)->ObtainReadData();
	
	char *transn = const_cast<char *> ("n");
  	char *transt = const_cast<char *> ("t");
    char *uplo = const_cast<char *> ("u");
    double one = 1, zero = 0, neg_one = -1;
    integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r; 
    
    double *UDptr = new double[MR];
    dgemm_(transn, transn, &M, &R, &R, &one, const_cast<double *> (Uptr), &M, const_cast<double *> (Dptr), &R, &zero, UDptr, &M);
	SharedSpace *Temp1 = new SharedSpace(2, m, n);
    double *Xptr = Temp1->ObtainWriteEntireData();
    dgemm_(transn, transt, &M, &N, &R, &one, UDptr, &M, const_cast<double *> (Vptr), &N, &zero, Xptr, &M);
	delete [] UDptr;
	    
    SharedSpace *Temp2 = new SharedSpace(2, m, n);
    double *Errptr = Temp2->ObtainWriteEntireData();
    dcopy_(&MN, A, &inc, Errptr, &inc);
    daxpy_(&MN, &neg_one, Xptr, &inc, Errptr, &inc);
	
	SharedSpace *Temp3 = new SharedSpace(2, m, n);
    double *QXptr = Temp3->ObtainWriteEntireData();
    dsymv_(uplo, &MN, &one, W, &MN, Errptr, &inc, &zero, QXptr, &inc);
	
	double result = 0;
    result = ddot_(&MN, Errptr, &inc, QXptr, &inc);
    if(UseGrad)
    {
        x->AddToTempData("X", Temp1);
        x->AddToTempData("err", Temp2);
        x->AddToTempData("QX", Temp3);
    }
    else
    {
        delete Temp1;
        delete Temp2;
        delete Temp3;
    }	
	return result;
};
//
//void WeightedLowRank::EucGrad(Variable *x, Vector *egf) const
//{
//	ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
//    double *Uptr = const_cast<double *> (ProdxxM->GetElement(0)->ObtainReadData());
//   	double *Dptr = const_cast<double *> (ProdxxM->GetElement(1)->ObtainReadData());
//    double *Vptr = const_cast<double *> (ProdxxM->GetElement(2)->ObtainReadData());
//	char *transn = const_cast<char *> ("n");
//  	char *transt = const_cast<char *> ("t");
//    double one = 1, zero = 0, neg_one = -1.0, neg_two = -2.0;
//    integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r, NR = n * r; 
//
//    const SharedSpace *Temp = x->ObtainReadTempData("QX");
//	const double *QXptr = Temp->ObtainReadData();
//	double *fullgrad = new double[MN];
//	dcopy_(&MN, const_cast<double *> (QXptr), &inc, fullgrad, &inc);
//	dscal_(&MN, &neg_two, fullgrad, &inc);
//
//	double *XiVptr = new double[MR];
//	dgemm_(transn, transn, &M, &R, &N, &one, fullgrad, &M, Vptr, &N, &zero, XiVptr, &M);
//	double *XiUptr = new double[NR];
//	dgemm_(transt, transn, &N, &R, &M, &one, fullgrad, &M, Uptr, &M, &zero, XiUptr, &N);
//	delete [] fullgrad;
//	
//	double *egfTV = egf->ObtainWriteEntireData();
//	double *Udotptr = egfTV;
//	double *Ddotptr = egfTV + m * r;
//	double *Vdotptr = Ddotptr + r * r;
//	dgemm_(transn, transt, &M, &R, &R, &one, XiVptr, &M, Dptr, &R, &zero, Udotptr, &M);
//	dgemm_(transt, transn, &R, &R, &M, &one, Uptr, &M, XiVptr, &M, &zero, Ddotptr, &R);
//	dgemm_(transn, transn, &N, &R, &R, &one, XiUptr, &N, Dptr, &R, &zero, Vdotptr, &N);
//	delete [] XiUptr;
//	delete [] XiVptr;
//};
//
//void WeightedLowRank::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
//{
//	etax->CopyTo(exix);
//};

void WeightedLowRank::RieGrad(Variable *x, Vector *gf) const        
{
    ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
    const double *Uptr = ProdxxM->GetElement(0)->ObtainReadData();
   	const double *Dptr = ProdxxM->GetElement(1)->ObtainReadData();
    const double *Vptr = ProdxxM->GetElement(2)->ObtainReadData();
	char *transn = const_cast<char *> ("n");
  	char *transt = const_cast<char *> ("t");
    double one = 1, zero = 0, neg_one = -1.0, neg_two = -2.0;
    integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r, NR = n * r, RR = r * r; 	
    const SharedSpace *Temp = x->ObtainReadTempData("QX");
	const double *QXptr = Temp->ObtainReadData();
	double *fullgrad = new double[MN];
	dcopy_(&MN, const_cast<double *> (QXptr), &inc, fullgrad, &inc);
	dscal_(&MN, &neg_two, fullgrad, &inc);
	//x->Print("x");//----
	//ForDebug::Print("fullgrad:", fullgrad, m, n);//----
	double *XiVptr = new double[MR];
	dgemm_(transn, transn, &M, &R, &N, &one, fullgrad, &M, const_cast<double *> (Vptr), &N, &zero, XiVptr, &M);
	
	double *XiUptr = new double[NR];
    dgemm_(transt, transn, &N, &R, &M, &one, fullgrad, &M, const_cast<double *> (Uptr), &M, &zero, XiUptr, &N);
	delete [] fullgrad;
	
	// compute inverse of D
	integer info;
	integer *IPIV = new integer[R+1];
	double *WORK = new double[RR];
	double *Dinv = new double[RR];
	dcopy_(&RR, const_cast<double *> (Dptr), &inc, Dinv, &inc);
	dgetrf_(&R, &R, Dinv, &R, IPIV, &info);
	dgetri_(&R, Dinv, &R, IPIV, WORK, &RR, &info);	
	delete [] IPIV;
	delete [] WORK;
	
	double *gfTV = gf->ObtainWriteEntireData();
	double *Udotptr = gfTV;
	double *Ddotptr = gfTV + m * r;
	double *Vdotptr = Ddotptr + r * r;
	dgemm_(transt, transn, &R, &R, &M, &one, const_cast<double *> (Uptr), &M, XiVptr, &M, &zero, Ddotptr, &R);
	dgemm_(transn, transn, &M, &R, &R, &one, const_cast<double *> (Uptr), &M, Ddotptr, &R, &zero, Udotptr, &M);
    dscal_(&MR, &neg_one, Udotptr, &inc);
    daxpy_(&MR, &one, XiVptr, &inc, Udotptr, &inc);
	dgemm_(transn, transt, &N, &R, &R, &one, const_cast<double *> (Vptr), &N, Ddotptr, &R, &zero, Vdotptr, &N);
    dscal_(&NR, &neg_one, Vdotptr, &inc);
    daxpy_(&NR, &one, XiUptr, &inc, Vdotptr, &inc);
	
	double *Udottemp = new double[MR];
	double *Vdottemp = new double[NR];
	dgemm_(transn, transn, &M, &R, &R, &one, Udotptr, &M, Dinv, &R, &zero, Udottemp, &M);
	dgemm_(transn, transt, &N, &R, &R, &one, Vdotptr, &N, Dinv, &R, &zero, Vdottemp, &N);
	dcopy_(&MR, Udottemp, &inc, Udotptr, &inc);
	dcopy_(&NR, Vdottemp, &inc, Vdotptr, &inc);
	//gf->Print("gf");//---
	delete [] Udottemp;
	delete [] Vdottemp;	
	delete [] Dinv;
	delete [] XiUptr;
	delete [] XiVptr;	
};

void WeightedLowRank::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
{
	etax->CopyTo(xix);	
	
	/*ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
	 double *Uptr = const_cast<double *> (ProdxxM->GetElement(0)->ObtainReadData());
	 double *Dptr = const_cast<double *> (ProdxxM->GetElement(1)->ObtainReadData());
	 double *Vptr = const_cast<double *> (ProdxxM->GetElement(2)->ObtainReadData());
	 char *transn = const_cast<char *> ("n");
	 char *transt = const_cast<char *> ("t");
	 double one = 1, zero = 0, neg_one = -1.0, neg_two = -2.0;
	 integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r, NR = n * r; 	
	 
	 const SharedSpace *Temp = x->ObtainReadTempData("QX");
	 const double *QXptr = Temp->ObtainReadData();
	 
	 double *gfTV = etax->ObtainWriteEntireData();
	 double *Udotptr = gfTV;
	 double *Ddotptr = gfTV + m * r;
	 double *Vdotptr = Ddotptr + r * r;
	 
	 double *UDdot = new double[MR];
	 dgemm_(transn, transn, &M, &R, &R, &one, Uptr, &M, Ddotptr, &R, &zero, UDdot, &M);
	 
	 double *Qeta = new double[MN];
	 dgemm_(transn, transt, &M, &N, &R, &one, Udotptr, &M, Vptr, &N, &zero, Qeta, &M);*/
	
	
	
	
	
	
};




