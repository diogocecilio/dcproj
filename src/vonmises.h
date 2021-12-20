#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"


class vonmises : elastoplasticbase
{
public:
	vonmises(Doub young,Doub nu, Doub sigy);
	vonmises();
	~vonmises();

	void closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, MatDoub & Dep, Doub & projgamma);
	Doub yield(Doub J2);
	MatDoub GetElasticMatrix();
	MatDoub GetInverseElasticMatrix();
	MatDoub F1HWCylVonMises(Doub xisol, Doub rho, Doub betasol);
	MatDoub HW(MatDoub sig);
	MatDoub stressrecosntruction(MatDoub val, MatDoub vec);
	MatDoub dadsig(NRtensor<Doub>  sigprojvoigt);
	MatDoub P();
	void SetElastic(){
         DebugStop();
    }
    void SetTangentMatrixType(bool type)
    {
         DebugStop();
    }
    Doub phi(NRtensor<Doub> epse)
    {
        NRmatrix<Doub>  tempepsemat, stresstrial;
        NRmatrix<Doub>  C = GetElasticMatrix();
        epse.FromTensorToNRmatrix(tempepsemat);
        C.Mult(tempepsemat, stresstrial);
        NRtensor<Doub>  stresstrialtensor;
        epse.FromNRmatrixToTensor(stresstrial, stresstrialtensor);
        Doub I1, J2;
        J2 = stresstrialtensor.J2();
        I1 = stresstrialtensor.I1();
        Doub xi = I1 / sqrt(3.);
        Doub rho = sqrt(2. * J2);
        Doub yieldcr = yield(J2);
        return yieldcr;
    }

	inline void setup(Doub young, Doub nu, Doub sigy) {
		fyoung = young;
		fsigy = sigy;
		fnu = nu;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
	}

	inline void GetElasticParameters(Doub young, Doub nu, Doub sigy, Doub K, Doub G)
	{
		young = fyoung;
		nu = fnu;
		sigy = fsigy;
		K = fK;
		G = fG;
	}
	inline void updateatributes(NRvector<MatDoub> mult)
	{
		fyoung0 = fyoung;
		fnu0 = fnu;
		fsigy0 = fsigy;
		Doub newyoung = fyoung + mult[0][0][0]*fyoung;
		setup(newyoung, fnu, fnu);
	}
	inline void restoreoriginalatributes()
	{
		setup(fyoung0, fnu0, fsigy0);
	}

	void GetMatConstants(NRvector<Doub>& consts) {
		consts.assign(3, 0.);
		consts[0] = fyoung;
		consts[1] = fnu;
		consts[2] = fsigy;
	}
	void SetMatConstants(NRvector<Doub>& consts) {
		setup(consts[0], consts[1], consts[2]);
	}

	inline void reset() {
		fyoung0 = 0.;
		fnu0 = 0.;
		fyoung = 0.;
		fnu = 0.;
		fG = 0.;
		fK = 0.;
	}

	Doub fyoung;
	Doub fnu;
	Doub fyoung0;
	Doub fnu0;
private:
	Doub fsigy0;
	Doub fsigy;
	Doub fK;
	Doub fG;
};

