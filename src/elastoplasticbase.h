#pragma once
#include "nr3.h"

class elastoplasticbase
{
public:

	~elastoplasticbase();

	virtual void closestpointproj(NRtensor<Doub> epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, MatDoub & Dep, Doub & projgamma)=0;
	//virtual Doub yield(Doub J2)=0;
	virtual MatDoub GetElasticMatrix()=0;
	virtual MatDoub GetInverseElasticMatrix()=0;
	//virtual MatDoub F1HWCylVonMises(Doub xisol, Doub rho, Doub betasol)=0;
	virtual MatDoub HW(MatDoub sig)=0;
	virtual MatDoub stressrecosntruction(MatDoub val, MatDoub vec)=0;
	virtual MatDoub dadsig(NRtensor<Doub>  sigprojvoigt)=0;
	virtual MatDoub P()=0;

	virtual void GetMatConstants(NRvector<Doub>& consts) = 0;
	virtual void SetMatConstants(NRvector<Doub>& consts) = 0;
//	virtual Doub GetCoes()=0;
//	virtual Doub GetPhi()=0;
//	virtual Doub GetYoung()=0;
//	virtual Doub GetNu()=0;
};

