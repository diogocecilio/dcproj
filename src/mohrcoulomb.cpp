#include "mohrcoulomb.h"

mohrcoulomb::mohrcoulomb()
{
    DebugStop();
}
mohrcoulomb::mohrcoulomb(Doub Phi, Doub Psi, Doub c,Doub young, Doub nu)
{
     DebugStop();
}
mohrcoulomb::mohrcoulomb(const mohrcoulomb &cp)
{
     DebugStop();
}

template<class T>
NRvector<T> mohrcoulomb::SigmaElastPV(const NRvector<T> &deform) const
{
    T trace = deform[0] + deform[1] + deform[2];
    NRvector<T> sigma(3, 0.);
    
    
    sigma = trace * flambda + 2 * fmu * deform[0];
    sigma = trace * flambda+ 2 *fmu * deform[1];
    sigma = trace * flambda + 2 * fmu * deform[2];
    
    return sigma;
}

template<class T>
T mohrcoulomb::PhiPlane(const NRvector<T> &sigma) const
{
    const Doub sinphi = sin(fPhi);
    const Doub cosphi = cos(fPhi);
    
    return sigma[0] - sigma[2] + (sigma[0] + sigma[2]) * sinphi - 2. * fc*cosphi;
}


template<class T>
bool mohrcoulomb::ReturnMapPlane(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const
{
    sigma_projected = sigma_trial;
    NRvector<T> eigenvalues = sigma_projected;
    const Doub sinphi = sin(fPhi);
    const Doub sinpsi = sin(fPsi);
    const Doub cosphi = cos(fPhi);
    const Doub sinphi2 = sinphi*sinphi;
    const Doub cosphi2 = 1. - sinphi2;
    const Doub constA = 4. * fmu *(1. + sinphi * sinpsi / 3.) + 4. * fK * sinphi*sinpsi;

    T phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * fc*cosphi;
    T gamma = memory.fGamma[0];
    int n_iterations = 30; 
    Doub H=0;
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {
        T jac = -constA - T(4. * cosphi2) * H;
        T delta_gamma = - phi / jac;
        gamma += delta_gamma;
        phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * fc * cosphi - constA * gamma;
        if (fabs(phi)<1.e-8) {
            break;
        }
    }

    if (i == n_iterations) {
        DebugStop();
    }
    
    Doub epsbar = H;
    memory.fGamma[0] = gamma;
    eigenvalues[0] -= T(2. * fmu*(1. + sinpsi / 3.) + 2. * fK * sinpsi) * gamma;
    eigenvalues[1] += T((4. * fmu / 3. - fK*2.) * sinpsi) * gamma;
    eigenvalues[2] += T(2. * fmu*(1 - sinpsi / 3.) - 2. * fK * sinpsi) * gamma;
    sigma_projected = eigenvalues;
    epsbarnew = epsbar;
    
    bool check_validity_Q = (eigenvalues[0] > eigenvalues[1]|| fabs(eigenvalues[0]-eigenvalues[1])<1.e-8 && (eigenvalues[1] > eigenvalues[2]) || (eigenvalues[1]-eigenvalues[2])<1.e-8);
    return (check_validity_Q);   
}

void mohrcoulomb::ComputePlaneTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    const Doub sin_phi = sin(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub G = fmu, K = fK;
    const Doub denominator = 6.0 * G + 2.0 * (G + 3.0 * K) * sin_phi * sin_psi;
    
    Doub epsbar = epsbarp;
    Doub c=fc, H=0;
    //PlasticityFunction(epsbar, c, H);
    
    tang.assign(3, 3,0.);
    
    // First column
    tang[0][0] = (sin_phi - 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    tang[1][0] = (2.0 * G - 3.0 * K) * (sin_phi + 1.0) * sin_psi / denominator;
    tang[2][0] = -(sin_phi + 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    
    // Second column
    tang[0][1] = 0.0;
    tang[1][1] = 1.0;
    tang[2][1] = 0.0;
    
    // Third column
    tang[0][2] = -(sin_phi - 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    tang[1][2] = (2.0 * G - 3.0 * K) * (sin_phi - 1.0) * sin_psi / denominator;
    tang[2][2] = (sin_phi + 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
}

template<class T>
bool mohrcoulomb::ReturnMapLeftEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const
{
    
}

void mohrcoulomb::ComputeLeftEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    
}

template<class T>
bool mohrcoulomb::ReturnMapRightEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const
{
    
}

void mohrcoulomb::ComputeRightEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    
}

template<class T>
bool mohrcoulomb::ReturnMapApex(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const
{
                
}
    

void mohrcoulomb::ComputeApexGradient(NRmatrix<Doub> & gradient, Doub & eps_bar_p) const
{
    
}
    

void mohrcoulomb::ProjectSigma(const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub> * gradient)
{
    
}


void mohrcoulomb::Phi(NRvector<Doub> sig_vec, Doub alpha, NRvector<Doub> &phi)const
{
    phi.resize(3);
    for (int i = 0; i < 3; i++) phi[i] = 0;
    phi[0] = PhiPlane(sig_vec);
    phi[2] = PhiPlane(sig_vec); // Consistency with two surfaces models
}
