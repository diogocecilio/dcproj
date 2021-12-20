#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"
#include "roots.h"
#include<Eigen/SparseCholesky>
#include <iostream>
#include <Eigen/Dense>

    /** @brief This class implements the Hyperbolic Drucker Prager elastoplastic model 
    Detailed description follows here.
    @author Diogo Cecílio
    @date dez 2021
    @brief Reference: <a href="https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29EM.1943-7889.0001737">link text</a>
    */

//using Eigen::MatrixXd;
using namespace Eigen;
class druckerprager
{
public:

    
    
    /**
     * @brief Default constructor
     */
	druckerprager();
    
    /**
	 * @brief Class constructor
	 * @param [in] young material youngs modulus
     * @param [in] nu material poisson ratio 
     * @param [in] coesion  material cohesion 
     * @param [in] frictionangle material internal friction angle
	 */
	druckerprager(Doub young, Doub nu, Doub coesion, Doub frictionangle);
    
    /**
     * @brief Default destructor
     */
	~druckerprager();

    /**
	 * @brief Implements the elastoplastic decomposition whit the Closest Point Projection Method and computes the tangent operator.
	 * @param [in] epst total strain tensor
     * @param [in] epsp plastic strain tensor
     * @param [out] projstress  projected stress tensor
     * @param [out] projstrain projected strain tensor
     * @param [out] Dep Elastoplastic tangent operator
     * @param [out] projgamma Plastic multiplier
	 */
	void closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma);
    
    /**
     * @brief Compute the yield function
     * @param [in] xi cylindrical coordinate xi=I1/3
     * @param [in] rho cylindrical coordinate rho =  sqrt(2 J2)
     * @param [out] yieldfunction 
     */
	Doub yield(Doub xi, Doub rho);
    
	NRmatrix<Doub> GetElasticMatrix();
	NRmatrix<Doub>  GetInverseElasticMatrix();
    
    /**
     * @brief Compute the cylindrical coordinates in DP plane strain match.
     * @param [in] xi cylindrical coordinate xi=I1/3
     * @param [in] rho cylindrical coordinate rho =  sqrt(2 J2)
     * @param [in] beta angle
     * @param [out] NRmatrix<Doub> {sig1,sig2,sig3} 
     */
	NRmatrix<Doub>  F1HWCylDruckerPragerSmoothPSMATCH(Doub xisol, Doub rho, Doub betasol);
    
    /**
     * @brief Compute the cylindrical coordinates
     * @param [in] sig principla stress
     * @param [out] NRmatrix<Doub> {xi,rho,beta}  cylindrical coordinates
     */
	NRmatrix<Doub>  HW(NRmatrix<Doub>  sig);
    
    /**
     * @brief Reconstruct a second order stress tensor
     * @param [in] val eigenvalues
     * @param [in] vec eigenvectors
     * @param [out] fulltensor 
     */
	NRmatrix<Doub>  stressrecosntruction(NRmatrix<Doub>  val, NRmatrix<Doub>  vec);
    
        
    /**
     * @brief Equation 58 of the Reference. Compute the flow rule derivative with respect to the projected stress 
     * @param [in] sigprojvoigt projected stress tensor in Voigt notation 
     * @param [out] dadsig derivative with respect to the projected stress 
     */
	NRmatrix<Doub>  dadsig(NRtensor<Doub>  sigprojvoigt);
    
    /**
     * @brief   \f$
     * \boldsymbol{P}=  \left[\begin{array}{cccccc}
     * 2/3 &  -1/3 & -1/3 & 0 & 0 & 0\\
     * -1/3 & 2/3 & -1/3 & 0 & 0 & 0\\
     * -1/3 & -1/3 & 2/3 & 0 & 0 & 0\\
     * 0 & 0 & 0 & 2 & 0 & 0\\
     * 0 & 0 & 0 & 0 & 2 & 0\\
     * 0 & 0 & 0 & 0 & 0 & 2\\
     * \end{array}\right]
     * \f$
     */
	NRmatrix<Doub>  P();
    /**
     * @brief Equation 57 of the Reference. Compute the flow rule
     * @param [in] sigprojvoigt projected stress tensor in Voigt notation 
     * @param [out] avec yield function derivative with respect to the projected stress 
     */
	NRmatrix<Doub>  avec(NRtensor<Doub>  sigprojvoigt);

    /**
     * @brief Find the project xi cylindrical coordinate using a combination of Newton-Raphson and bisection.
     * @param [in] pt trial principal stress  
     * @param [in] xitrial trial xi cylindrical coordinate
     */
	Doub FindMinimum(NRmatrix<Doub>  pt,Doub xitrial,bool flag);
    
    /**
     * @brief Compute  \f$ \boldsymbol{Q} =  \boldsymbol{I} + \Delta \gamma \mathbf{D}^e \frac{\boldsymbol{\partial{a}}}{\boldsymbol{\partial{\sigma}}}  \f$
     * @param [in] pt trial principal stress  
     * @param [in] xitrial trial xi cylindrical coordinate
     */
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
	
	//true consistent tangent Dep
	//false pseudo tangent R
	void SetTangentMatrixType(bool type)
    {
         fsetconsistentangent = type;
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
	bool fflag = true;//Initialised as  consistent tangent Dep
	bool fsetconsistentangent=true;
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


