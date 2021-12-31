#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include <iostream>

#include <random>


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
    Doub ftol=1.e-8;
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

    NRmatrix<Doub>  stressrecosntruction(NRvector<Doub>  val, NRmatrix<Doub>  vec)
    {
        NRmatrix<Doub>  sol(3, 3, 0.);
        for (Int i = 0;i < 3;i++)
        {
            MatDoub colvec(3, 1, 0.), colvect, temp;
            for (Int j = 0;j < 3;j++)
            {
                colvec[j][0] = vec[j][i];
            }

            colvec.Transpose(colvect);
            colvec.Mult(colvect, temp);
            temp *= val[i];
            sol += temp;
        }
        //sol.Print();
        return sol;
    }
    
    NRmatrix<Doub>  HW(NRmatrix<Doub>  sig)
{
	MatDoub sol(3, 1);
	Doub xi, rho, beta;
	xi = sig[0][0];
	rho = sig[1][0];
	beta = sig[2][0];
	sol[0][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta);
	sol[1][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta - 2. *  M_PI / 3.);
	sol[2][0] = xi / sqrt(3.) + sqrt(2. / 3.)*rho*cos(beta + 2. *  M_PI / 3.);
	return sol;
}
    
    
    void SetUp(Doub Phi, Doub Psi, Doub c,Doub young, Doub nu) {
        fPhi = Phi;
        fPsi = Psi;
        fc = c;
        fyoung = young;
        fnu = nu;
        fmu = fyoung/(2. * (1. + fnu));
        flambda = fyoung*fnu/((1. + fnu)* (1. - 2.* fnu));
        //cout << "YOUNG = " << fyoung << endl;
        //cout << "fnu" << fnu << endl;
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

       bool  ProjectSigma(NRvector<Doub> &sigma_trial, Doub &k_prev,NRvector<Doub> &sigma,Doub &k_proj);
    
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

			void ComputePlaneTangent(NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order) const;

    template<class T>
    bool ReturnMapLeftEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;

			void ComputeLeftEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order) const;

    template<class T>
    bool ReturnMapRightEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;

			void ComputeRightEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order) const;

    template<class T>
    bool ReturnMapApex(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const;
    

			void ComputeApexGradient(NRmatrix<Doub> & gradient, Doub & eps_bar_p, NRvector<Int> order) const;
    

			void ProjectSigma(const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub>  &gradient,NRvector<Int> order ,NRmatrix<Doub> &eigenvec,NRvector<Doub> pstrain);


    void Phi(NRvector<Doub> sig_vec, Doub alpha, NRvector<Doub> &phi)const;

    Doub Phi() {
        return fPhi;
    }
    
            
                void ComputePlanePrincipalStressDeriv(NRmatrix<Doub> &DPSTRS, Doub &epsbarp,NRvector<Int> order) const;
	
				void ComputePlaneDep(NRvector<Doub> strain_trial,NRvector<Doub> sigma_proj,NRmatrix<Doub> &eigenvec, Doub &epsbarp,NRvector<Int> order,NRmatrix<Doub> &Dep) const;
	
	
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
        Doub alpha=0;
        
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
	void closestpointproj2(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma);
	
	/**
	 * @brief Implements the elastoplastic decomposition whit the Closest Point Projection Method and computes the tangent operator.
	 * @param [in] epst total strain tensor
	 * @param [in] epsp plastic strain tensor
	 * @param [out] projstress  projected stress tensor
	 * @param [out] projstrain projected strain tensorr
	 */
	bool closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain);
	
	void closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp,NRtensor<Doub> &projstress, NRtensor<Doub> &projstrain,NRmatrix<Doub> &Dep2, Doub & projgamma)
	{
		projgamma=0;
		Dep2.assign(3,3,0.);
		NRmatrix<Doub>Dep(6,6,0.);
		srand((unsigned int)time(NULL));
		float a = 1.e-8;
		for(Int i=0;i<6;i++)
		{
			//if(fabs(epst[i])<1.e-20)epst[i]=1.e-20;
		}
		
		NRtensor<Doub> epstpertub(epst);
		
		bool iselastic = closestpointproj(epst, epsp, projstress, projstrain);
		if(iselastic==false)
		{
			NRvector<NRtensor<Doub>> vecprojperturbstress(6),vecprojperturbstrain(6);
			for(Int iperturb=0;iperturb<6;iperturb++)
			{
				NRtensor<Doub> epstpertub(epst), projstressperturb, projstrainperturb;
				epstpertub[iperturb]+=epst[iperturb]*a;
				//closestpointproj(epst, epsp, projstress, projstrain);
				closestpointproj(epstpertub, epsp, projstressperturb, projstrainperturb);

				NRtensor<Doub> dsig(projstressperturb),deps(epstpertub);
			
				dsig-=projstress;
				deps-=epst;
			
				for(Int i=0;i<6;i++)
				{
					if(fabs(deps[i])<1.e-30)deps[i]=1.e-30;
					if(fabs(dsig[i])<1.e-30)dsig[i]=1.e-30;
				}
			
				for(Int ivar=0;ivar<6;ivar++)
				{
					Dep[ivar][iperturb]=dsig[ivar]/deps[iperturb];
				//Dep[ivar][iperturb]=dsig[ivar]/a;
				}
			
			}
		}else{
				Dep = GetElasticMatrix();
			}
			
		Dep2[0][0] = Dep[0][0];Dep2[0][1] = Dep[0][1];Dep2[0][2] = Dep[0][5];
		Dep2[1][0] = Dep[1][0];Dep2[1][1] = Dep[1][1];Dep2[1][2] = Dep[1][5];
		Dep2[2][0] = Dep[5][0];Dep2[2][1] = Dep[5][1];Dep2[2][2] = Dep[5][5];
	
	}
	
	
	
	void ComputeNumericalDep2(NRtensor<Doub>  epst, NRtensor<Doub>  epsp,NRtensor<Doub> &projstress, NRtensor<Doub> &projstrain,NRmatrix<Doub> &Dep)
	{
		
		srand((unsigned int)time(NULL));
		float a = 1.e-8;
 
		epst.Print();
		
		NRtensor<Doub> randtensor,projstressperturb,  projstrainperturb;
		
		randtensor.XX()=(float(rand())/float((RAND_MAX)) );
		randtensor.YY()=(float(rand())/float((RAND_MAX)) );
		randtensor.ZZ()=(float(rand())/float((RAND_MAX)) );
		randtensor.XZ()=(float(rand())/float((RAND_MAX)) );
		randtensor.YZ()=(float(rand())/float((RAND_MAX)) );
		randtensor.XY()=(float(rand())/float((RAND_MAX)) );
		NRtensor<Doub> epstpertub(epst);
		//for(Int i=0;i<6;i++)epstpertub[i]=epst[i]+randtensor[i] *a;
		//for(Int i=0;i<6;i++)epstpertub[i]=epst[i]+epst[i]*a;
		//epstpertub+=a;
		epstpertub.XY()+=epst.XY()*a;
		//epstpertub.Print();
		//dsigproj/depstrial
		//dsigxx/depstxx dsigxx/depstyy dsigxx/depstzz dsigxx/depstxz dsigxx/depstyz dsigxx/depstxy
		//dsigyy/depstxx dsigyy/depstyy dsigyy/depstzz dsigyy/depstxz dsigyy/depstyz dsigyy/depstxy
		//dsigzz/depstxx dsigzz/depstyy dsigzz/depstzz dsigzz/depstxz dsigzz/depstyz dsigzz/depstxy
		//dsigxz/depstxx dsigxz/depstyy dsigxz/depstzz dsigxz/depstxz dsigxz/depstyz dsigxz/depstxy
		//dsigyz/depstxx dsigyz/depstyy dsigyz/depstzz dsigyz/depstxz dsigyz/depstyz dsigyz/depstxy
		//dsigxy/depstxx dsigxy/depstyy dsigxy/depstzz dsigxy/depstxz dsigxy/depstyz dsigxy/depstxy
		
		// dsig/deps = (sig(eps+ alpha d)-sig(eps))/alpha para alpha->0
		//alpha =1.e-6
		//d = random direction
		// dsigxx/depstxx = (sigpertubxx-sigxx)/(epstpertubxx-epstxx) ...
		
		closestpointproj(epst, epsp, projstress, projstrain);
		cout << "projstress"<<endl;
		projstress.Print();
		closestpointproj(epstpertub, epsp, projstressperturb, projstrainperturb);
		
		Dep.assign(6,6,0.);
		NRtensor<Doub> dsig(projstressperturb),deps(epstpertub);
		cout << "projstressperturb"<<endl;
		projstressperturb.Print();
		
		
		dsig-=projstress;
		deps-=epst;
		
		//for(Int i=0;i<6;i++)deps[i]=a;
Dep[0][0]=dsig.XX()/deps.XX();Dep[0][1]=dsig.XX()/deps.YY(); Dep[0][2]=dsig.XX()/deps.ZZ(); Dep[0][3]=dsig.XX()/deps.XZ(); Dep[0][4]=dsig.XX()/deps.YZ(); Dep[0][5]=dsig.XX()/deps.XY();

Dep[1][0]=dsig.YY()/deps.XX();Dep[1][1]=dsig.YY()/deps.YY(); Dep[1][2]=dsig.YY()/deps.ZZ(); Dep[1][3]=dsig.YY()/deps.XZ(); Dep[1][4]=dsig.YY()/deps.YZ(); Dep[1][5]=dsig.YY()/deps.XY();

Dep[2][0]=dsig.ZZ()/deps.XX();Dep[2][1]=dsig.ZZ()/deps.YY(); Dep[2][2]=dsig.ZZ()/deps.ZZ(); Dep[2][3]=dsig.ZZ()/deps.XZ(); Dep[2][4]=dsig.ZZ()/deps.YZ(); Dep[2][5]=dsig.ZZ()/deps.XY();

Dep[3][0]=dsig.XZ()/deps.XX();Dep[3][1]=dsig.XZ()/deps.YY(); Dep[3][2]=dsig.XZ()/deps.ZZ(); Dep[3][3]=dsig.XZ()/deps.XZ(); Dep[3][4]=dsig.XZ()/deps.YZ(); Dep[3][5]=dsig.XZ()/deps.XY();

Dep[4][0]=dsig.YZ()/deps.XX();Dep[4][1]=dsig.YZ()/deps.YY(); Dep[4][2]=dsig.YZ()/deps.ZZ(); Dep[4][3]=dsig.YZ()/deps.XZ(); Dep[4][4]=dsig.YZ()/deps.YZ(); Dep[4][5]=dsig.YZ()/deps.XY();

Dep[5][0]=dsig.XY()/deps.XX();Dep[5][1]=dsig.XY()/deps.YY(); Dep[5][2]=dsig.XY()/deps.ZZ(); Dep[5][3]=dsig.XY()/deps.XZ(); Dep[5][4]=dsig.XY()/deps.YZ(); Dep[5][5]=dsig.XY()/deps.XY();
	}
    
    
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


