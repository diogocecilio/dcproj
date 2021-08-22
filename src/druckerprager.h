#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"
#include "roots.h"
#include<Eigen/SparseCholesky>
#include <iostream>
#include <Eigen/Dense>

//using Eigen::MatrixXd;
using namespace Eigen;
class druckerprager
{
public:
	druckerprager(Doub young, Doub nu, Doub coesion, Doub frictionangle);
	druckerprager();
	~druckerprager();

	void closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma);
	Doub yield(Doub xi, Doub rho);
	NRmatrix<Doub> GetElasticMatrix();
	NRmatrix<Doub>  GetInverseElasticMatrix();
	NRmatrix<Doub>  F1HWCylDruckerPragerSmoothPSMATCH(Doub xisol, Doub rho, Doub betasol);
	NRmatrix<Doub>  HW(NRmatrix<Doub>  sig);
	NRmatrix<Doub>  stressrecosntruction(NRmatrix<Doub>  val, NRmatrix<Doub>  vec);
	NRmatrix<Doub>  dadsig(NRtensor<Doub>  sigprojvoigt);
	NRmatrix<Doub>  P();
	NRmatrix<Doub>  avec(NRtensor<Doub>  sigprojvoigt);

	Doub FindMinimum(NRmatrix<Doub>  pt,Doub xitrial,bool flag);
	MatDoub ComputeQ(NRmatrix<Doub>  fulltensorproj, NRmatrix<Doub>  tempepsemat, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, Doub & projgamma,NRmatrix<Doub>  & nvec);

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
        Doub yieldcr =-(pow(((fapex - xi / sqrt(3.)) / fa), 2) - pow(((rho / sqrt(2.)) / fb), 2) - 1);
        return yieldcr;
    }
	
	// Doub func2(Doub xi)

	//{
	//	Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
	//	Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

	//	Doub fun = ((4. * pow(sqrt(3.)*sig1 + sqrt(3)*sig2 + sqrt(3.)*sig3 - 3. * xi, 2)) / fK + (9. * pow(-2. * sig1 + sig2 + sig3 +
	//		2 * sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2.))*cos(beta), 2.)) / fG +
	//		(3 * pow(-3. * sig2 + 3. * sig3 + 2 * sqrt(3.)*sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2))*sin(beta), 2)) / fG) / 108.;
	//	return fun;
	//}
	//Doub   distfunddp(MatDoub &pt, VecDoub_I &x);
	inline void setup(Doub young, Doub nu, Doub coesion, Doub frictionangle) {
		fyoung = young;
		fnu = nu;
		fcoesion = coesion;
		fphi = frictionangle;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
		ftanphi = (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		fapex = fcoesion * 1. / tan(fphi);
		fa = (fcoesion / (sqrt(3)*tan(fphi)) - fapex);
		fb = fa*ftanphi;

	}

	Doub GetCoes()
	{
		return fcoesion;
	}

	Doub GetPhi()
	{
		return fphi;
	}

	Doub GetYoung()
	{
		return fyoung;
	}

	Doub GetNu()
	{
		return fnu;
	}

	inline void reset() {
		fyoung0 = 0.;
		fnu0 = 0.;
		fcoesion0 = 0.;
		fphi0 = 0.;
		fyoung = 0.;
		fnu = 0.;
		fcoesion = 0.;
		fphi = 0.;
		fG = 0.;
		fK = 0.;
		ftanphi = 0.;
		fapex = 0.;
		fa =0.;
		fb = 0.;
	}

	void updateatributes(NRvector<MatDoub> mult);

	inline void restoreoriginalatributes()
	{
		setup(fyoung0, fnu0, fcoesion0, fphi0);
	}


	void GetMatConstants(NRvector<Doub>& consts) {
		consts.assign(4, 0.);
		consts[0] = fyoung;
		consts[1] = fnu;
		consts[2] = fcoesion;
		consts[3] = fphi;
		
	}

	void SetMatConstants(NRvector<Doub>& consts) {
		setup(consts[0], consts[1], consts[2], consts[3]);
	}



	inline Doub Det(NRmatrix<Doub>  A)
	{

		MatrixXd AA(A.nrows(), A.nrows());
		for (int i = 0;i < A.nrows();i++)
		{
			for (int j = 0;j < A.ncols();j++)
			{
				AA(i, j) = A[i][j];
			}
		}
		Doub det = AA.determinant();
		return det;
	}

	inline NRmatrix<Doub>  Inverse(NRmatrix<Doub>  A)
	{
		MatDoub InvOut(A.nrows(), A.nrows());
		MatrixXd AA(A.nrows(), A.nrows()),Inv;
		for (int i = 0;i < A.nrows();i++)
		{
			for (int j = 0;j < A.ncols();j++)
			{
				AA(i, j) = A[i][j];
			}
		}
		Inv = AA.inverse();
		for (int i = 0;i < A.nrows();i++)
		{
			for (int j = 0;j < A.ncols();j++)
			{
				InvOut[i][j]= Inv(i, j);
			}
		}
		return InvOut;
	}


	Doub fyoung;
	Doub fnu;

private:
	Doub fcoesion;
	Doub fphi;
	Doub fcoesion0;
	Doub fphi0;

	Doub fK;
	Doub fG;
	Doub ftanphi;
	Doub fapex;
	Doub fa, fb;
	bool fflag = true;
	Doub fyoung0;

	Doub fnu0;
};


struct Funcstruct:druckerprager {

	NRmatrix<Doub>  pt;
	Doub fyoung;
	Doub fnu;
	Doub fcoesion;
	Doub fphi;
	Doub fK;
	Doub fG;
	Doub ftanphi;
	Doub fapex;
	Doub fa, fb;



	Doub operator()(VecDoub_I &x)
	{
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x[0];
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));

		Doub func2 = ((4. * pow(sqrt(3.)*sig1 + sqrt(3)*sig2 + sqrt(3.)*sig3 - 3. * xi, 2)) / fK + (9. * pow(-2. * sig1 + sig2 + sig3 +
			2 * sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2.))*cos(beta), 2.)) / fG +
			(3 * pow(-3. * sig2 + 3. * sig3 + 2 * sqrt(3.)*sqrt((pow(fb, 2)*(-3. * pow(fa, 2) + 3. * pow(fapex, 2) - 2. * sqrt(3.)*fapex*xi + pow(xi, 2))) / pow(fa, 2))*sin(beta), 2)) / fG) / 108.;
		return func2;
	}

	Doub operator()(const Doub &x)
	{
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x;
		if (xi > fapex) {
		//	xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG, K = fK, a = fa, b = fb, apex = fapex;
		Doub func2 = ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
		return func2;
	}

	void df(VecDoub_I &x, VecDoub_O &deriv)
	{
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x[0];
		if (xi > fapex) {
		//	xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG,K = fK, a = fa, b = fb, apex = fapex;
		deriv[0] = ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}


	void setpt(NRmatrix<Doub>  ptext, Doub young, Doub nu, Doub coesion, Doub frictionangle)
	{
		pt = ptext;
		setup(young,nu,coesion,frictionangle);
	}
	inline void setup(Doub young, Doub nu, Doub coesion, Doub frictionangle) {
		fyoung = young;
		fnu = nu;
		fcoesion = coesion;
		fphi = frictionangle;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
		ftanphi = (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		fapex = fcoesion * 1. / tan(fphi);


		//fa= (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		//fb = (3.) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		//fa = (6.*sin(fphi)) / ( sqrt(3.)*(3 + sin(fphi)) );
		//fb = (6.*cos(fphi)) / (sqrt(3.)*(3 + sin(fphi)));

		//fa = (6.*sin(fphi)) / (sqrt(3.)*(3 - sin(fphi)));//1.10
		//fb = (6.*cos(fphi)) / (sqrt(3.)*(3 - sin(fphi)));

		fa = (fcoesion / (sqrt(3)*tan(fphi)) - fapex);
		fb = fa*ftanphi;

	}
};


struct Funcd :druckerprager {
	MatDoub pt;
	Doub fyoung;
	Doub fnu;
	Doub fcoesion;
	Doub fphi;
	Doub fK;
	Doub fG;
	Doub ftanphi;
	Doub fapex;
	Doub fa, fb;
	Doub operator() (const Doub x) {
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x;
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG, K = fK, a = fa, b = fb, apex = fapex;
		return ((-24 * (sqrt(3)*sig1 + sqrt(3)*sig2 + sqrt(3)*sig3 - 3 * xi)) / K + (18 * pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*cos(beta)*
			(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(6 * sqrt(3)*pow(b, 2)*(-2 * sqrt(3)*apex + 2 * xi)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}
	Doub df(const Doub x) {
		Doub sig1 = pt[0][0], sig2 = pt[0][1], sig3 = pt[0][2];
		Doub xi = x;
		if (xi > fapex) {
			//xi = fapex;
		}
		Doub beta = atan((sqrt(3)*(-sig2 + sig3)) / (-2 * sig1 + sig2 + sig3));
		Doub G = fG, K = fK, a = fa, b = fb, apex = fapex;
		return (72 / K + (72 * pow(b, 2)*pow(-(sqrt(3)*apex) + xi, 2)*pow(cos(beta), 2)) / (pow(a, 2)*G*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) -
			(9 * pow(b, 4)*pow(2 * sqrt(3)*apex - 2 * xi, 2)*cos(beta)*(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 4)*G*pow((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2), 1.5)) +
			(36 * pow(b, 2)*cos(beta)*(-2 * sig1 + sig2 + sig3 + 2 * sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*cos(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))) +
			(72 * pow(b, 2)*pow(-(sqrt(3)*apex) + xi, 2)*pow(sin(beta), 2)) / (pow(a, 2)*G*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) -
			(3 * sqrt(3)*pow(b, 4)*pow(2 * sqrt(3)*apex - 2 * xi, 2)*sin(beta)*(-3 * sig2 + 3 * sig3 +
				2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 4)*G*pow((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2), 1.5)) +
			(12 * sqrt(3)*pow(b, 2)*sin(beta)*(-3 * sig2 + 3 * sig3 + 2 * sqrt(3)*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2))*sin(beta))) /
			(pow(a, 2)*G*sqrt((pow(b, 2)*(-3 * pow(a, 2) + 3 * pow(apex, 2) - 2 * sqrt(3)*apex*xi + pow(xi, 2))) / pow(a, 2)))) / 108.;
	}

	void setpt(NRmatrix<Doub>  ptext, Doub young, Doub nu, Doub coesion, Doub frictionangle)
	{
		pt = ptext;
		setup(young, nu, coesion, frictionangle);
	}
	inline void setup(Doub young, Doub nu, Doub coesion, Doub frictionangle) {
		fyoung = young;
		fnu = nu;
		fcoesion = coesion;
		fphi = frictionangle;
		fG = young / (2.* (1. + nu));
		fK = young / (3.* (1. - 2.* nu));
		ftanphi = (3.*tan(fphi)) / sqrt(9. + 12.* tan(fphi) *tan(fphi));
		fapex = fcoesion * 1. / tan(fphi);
		fa = (fcoesion / (sqrt(3)*tan(fphi)) - fapex);
		fb = fa*ftanphi;

	}

};


