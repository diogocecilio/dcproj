#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include <iostream>




class mohrcoulomb  {  
public:

    enum {
        NYield = 3
    };

private:
    
    Doub fPhi;
    Doub fPsi;
    Doub fc;
    Doub fyoung;
    Doub fnu;
    Doub flambda;
    Doub fmu;
    Doub fK;
    
    Doub fyoung0 ;
	Doub fnu0 ;
	Doub fc0;
	Doub fPhi0;
    Doub fPsi0 ;
    Doub ftol=1.e-5;
    bool fsetconsistentangent=true;


protected:
    Doub fEpsPlasticBar;

public:

    /// structure which contains the decision tree of the return map
    // we can only expect a consistent tangent matrix if the decision tree remains the same

    struct TComputeSequence {

        TComputeSequence() : fWhichPlane(ENoPlane), fGamma(0) {

        }

        TComputeSequence(const TComputeSequence &copy) : fWhichPlane(copy.fWhichPlane), fGamma(copy.fGamma) {

        }

        TComputeSequence &operator=(const TComputeSequence &copy) {
            fWhichPlane = copy.fWhichPlane;
            fGamma = copy.fGamma;
            return *this;
        }

        enum MPlane {
            ENoPlane, EElastic, EMainPlane, ERightEdge, ELeftEdge, EApex
        };

        MPlane fWhichPlane;

        NRvector<Doub> fGamma;
    };

public:

    mohrcoulomb();
    mohrcoulomb(Doub Phi, Doub Psi, Doub c,Doub young, Doub nu);
    mohrcoulomb(const mohrcoulomb &cp);

    
    
    
    void SetUp(Doub Phi, Doub Psi, Doub c,Doub young, Doub nu) {
        fPhi = Phi;
        fPsi = Psi;
        fc = c;
        fyoung = young;
        fnu = nu;
        fmu = fyoung/(2. * (1. + fnu));
        flambda = fyoung*fnu/((1. + fnu)* (1. - 2.* fnu));
        fK = fyoung / (3.* (1. - 2.* fnu));
    }

    mohrcoulomb & operator=(const mohrcoulomb &cp);

    void SetEpsBar(Doub &epsbar) {
        fEpsPlasticBar = epsbar;
    }


    void SetElasticResponse(Doub young, Doub nu) {
        fyoung = young;
        fnu = nu;
    }

    
    Doub  InitialDamage(const NRvector<Doub> &stress_p)const;
    
    template <class T>
    void PlasticityFunction(const T epsp, T &c, T &H) const;
    
    template<class T>
    NRvector<T> SigmaElastPV(const NRvector<T> &deform) const;

    template<class T>
    T PhiPlane(const NRvector<T> &sigma) const;


    template<class T>
    bool ReturnMapPlane(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;

    void ComputePlaneTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const;

    template<class T>
    bool ReturnMapLeftEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;

    void ComputeLeftEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const;

    template<class T>
    bool ReturnMapRightEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;

    void ComputeRightEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const;

    template<class T>
    bool ReturnMapApex(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;
    

    void ComputeApexGradient(NRmatrix<Doub> & gradient, Doub & eps_bar_p) const;
    

    void ProjectSigma(const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub>  &gradient );


    void Phi(NRvector<Doub> sig_vec, Doub alpha, NRvector<Doub> &phi)const;

    Doub Phi() {
        return fPhi;
    }
    
    
    NRvector<Doub>  phi(NRtensor<Doub> epse)
    {
        NRmatrix<Doub>  tempepsemat, stresstrial,pt,vec;
        NRmatrix<Doub>  C = GetElasticMatrix();
        epse.FromTensorToNRmatrix(tempepsemat);
        C.Mult(tempepsemat, stresstrial);
        NRtensor<Doub>  stresstrialtensor;
        epse.FromNRmatrixToTensor(stresstrial, stresstrialtensor);
        
        NRvector<Doub>  sig_vec,phiv ;
 
        stresstrialtensor.EigenSystem(pt, vec);
    
        sig_vec[0] = pt[0][0];
        sig_vec[1] =  pt[0][1];
        sig_vec[2] = pt[0][2];
        Doub alpha;
        
        Phi(sig_vec, alpha, phiv);
        
        return phiv;
    }
    
    void updateatributes(NRvector<MatDoub> mult);
    
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
    
    
    /// Set up the phi
    void SetPhi(Doub phi)
    {
        fPhi = phi;
    }
    
    Doub Psi() {
        return fPsi;
    }


    Doub Cohesion() {
        return fc;
    }
    
    void SetCohesion(Doub cohesion)
    {
        fc = cohesion;
    }

    
    Doub E() {
        return fyoung;
    }

    Doub Poisson() {
        return fnu;
    }
    
    void GetMatConstants(NRvector<Doub>& consts) {
		consts.assign(5, 0.);
		consts[0] = fyoung;
		consts[1] = fnu;
		consts[2] = fc ;
		consts[3] = fPhi;
        consts[4] = fPsi;
		
	}
	
    inline void restoreoriginalatributes()
	{
		SetUp(fPhi0, fPsi0, fc0, fyoung0,fnu0);
	}
	
		void SetMatConstants(NRvector<Doub>& consts) {
		SetUp(consts[0], consts[1], consts[2], consts[3],consts[4]);
	}

		void SetTangentMatrixType(bool type)
    {
         fsetconsistentangent = type;
    }
NRmatrix<Doub>  GetElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fmu, K = fK;
	C[0][0] = (4 * G) / 3 + K;    C[0][1] = -((2 * G) / 3) + K;C[0][2] = -((2 * G) / 3) + K;C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -((2 * G) / 3) + K; C[1][1] = (4 * G) / 3 + K;C[1][2] = -((2 * G) / 3) + K;C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -((2 * G) / 3) + K; C[2][1] = -((2 * G) / 3) + K;C[2][2] = (4 * G) / 3 + K;C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                  C[3][1] = 0;                 C[3][2] = 0;                 C[3][3] = G; C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                  C[4][1] = 0;                 C[4][2] = 0;                 C[4][3] = 0.;C[4][4] = G; C[4][5] = 0.;
	C[5][0] = 0;                  C[5][1] = 0;                 C[5][2] = 0;                 C[5][3] = 0.;C[5][4] = 0.;C[5][5] = G;
	return C;
}
NRmatrix<Doub>   GetInverseElasticMatrix()
{
	MatDoub C(6, 6, 0.);
	Doub G = fmu, K = fK;
	C[0][0] = (G + 3 * K) / (9.*G*K);    C[0][1] = -1 / (6.*G) + 1 / (9.*K);C[0][2] = -1 / (6.*G) + 1 / (9.*K);C[0][3] = 0.;C[0][4] = 0.;C[0][5] = 0.;
	C[1][0] = -1 / (6.*G) + 1 / (9.*K);  C[1][1] = (G + 3 * K) / (9.*G*K);  C[1][2] = -1 / (6.*G) + 1 / (9.*K);C[1][3] = 0.;C[1][4] = 0.;C[1][5] = 0.;
	C[2][0] = -1 / (6.*G) + 1 / (9.*K);  C[2][1] = -1 / (6.*G) + 1 / (9.*K);C[2][2] = (G + 3 * K) / (9.*G*K);  C[2][3] = 0.;C[2][4] = 0.;C[2][5] = 0.;
	C[3][0] = 0;                         C[3][1] = 0;                       C[3][2] = 0;                       C[3][3] = 1. / G;C[3][4] = 0.;C[3][5] = 0.;
	C[4][0] = 0;                         C[4][1] = 0;                       C[4][2] = 0;                       C[4][3] = 0.;C[4][4] = 1. / G;C[4][5] = 0.;
	C[5][0] = 0;                         C[5][1] = 0;                       C[5][2] = 0;                       C[5][3] = 0.;C[5][4] = 0.;C[5][5] = 1. / G;
	return C;
}

};


