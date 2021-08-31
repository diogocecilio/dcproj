#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"
#include "roots.h"
#include<Eigen/SparseCholesky>
#include <iostream>
#include <Eigen/Dense>



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
    

    void ProjectSigma(const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub> * gradient = NULL);


    void Phi(NRvector<Doub> sig_vec, Doub alpha, NRvector<Doub> &phi)const;

    Doub Phi() {
        return fPhi;
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

};


