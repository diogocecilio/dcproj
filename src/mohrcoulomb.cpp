#include "mohrcoulomb.h"

mohrcoulomb::mohrcoulomb()
{
    DebugStop();
}
mohrcoulomb::mohrcoulomb(Doub Phi, Doub Psi, Doub c,Doub young, Doub nu)
{
     fPhi = Phi;
     fPsi = Psi;
     fc = c;
     fyoung = young;
     fnu = nu;
     SetUp( Phi,  Psi,  c, young,  nu);
}
mohrcoulomb::mohrcoulomb(const mohrcoulomb &cp)
{
     DebugStop();
}

Doub mohrcoulomb::InitialDamage(const NRvector<Doub> &stress_p) const{
    
    std::cout << "TPZYCMohrCoulombPV::There is no damage variable for this model at the current time." << std::endl;
    
    NRvector<Doub>  phi(3);
    Doub alpha = 0.0;
    Doub tol = 1.0e-10;
    Phi(stress_p, alpha, phi);
    bool Is_valid_stress_on_cap_Q =  fabs(phi[0]) < tol || phi[0] < 0.0;
    
    if (!Is_valid_stress_on_cap_Q) {
        std::cerr << "MohrCoulombPV::Invalid stress state." << std::endl;
        DebugStop();
    }

    return 0.0;
}

void mohrcoulomb::Phi(NRvector<Doub>sig_vec, Doub alpha, NRvector<Doub> &phi)const {
    phi.resize(3);
    for (int i = 0; i < 3; i++) phi[i] = 0;
    phi[0] = PhiPlane(sig_vec);
    phi[2] = PhiPlane(sig_vec); // Consistency with two surfaces models
}

template <class T>
void mohrcoulomb::PlasticityFunction(const T epsp, T &c, T &H) const {
    c = fc; // c(epsp)
    H = 0.; // dc(epsp)/depsp
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
    Doub temp =sigma[0] - sigma[2] + (sigma[0] + sigma[2]) * sinphi;
    cout << "sigma[0]" << sigma[0] <<endl;
    cout << "sigma[2]" << sigma[2] <<endl;
    cout << "temp" << temp <<endl;
    return  temp - 2. * fc*cosphi;
}
//YOUNG =    20000.000000000000       POISS =   0.48999999999999999       SINPHI =   0.34202014332566871     
// PSTRS1 =    8548.2496567310754     
// PSTRS2 =    7820.6523120359971     
// PSTRS3 =    7412.2652657913695     
// SMCT =    6594.8019922923067     
// COHE =    50.000000000000000     
// R2CPHI =    1.8793852415718169     
// PHIA =    6500.8327302137159     
// MAIN PLANE
// S1 =   -42.313073669229198       S2 =   -133.00065302404710       S3 =   -229.11683046147573

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
    T c, H;
    T epsbar = T(fEpsPlasticBar + memory.fGamma[0]*2. * cosphi);
    PlasticityFunction(epsbar, c, H);
    
    T phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * fc*cosphi;
    T gamma = memory.fGamma[0];
    int n_iterations = 30; 

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
    
    //Doub epsbar = H;
    memory.fGamma[0] = gamma;
    eigenvalues[0] -= T(2. * fmu*(1. + sinpsi / 3.) + 2. * fK * sinpsi) * gamma;
    eigenvalues[1] += T((4. * fmu / 3. - fK*2.) * sinpsi) * gamma;
    eigenvalues[2] += T(2. * fmu*(1 - sinpsi / 3.) - 2. * fK * sinpsi) * gamma;
    sigma_projected = eigenvalues;
    sigma_projected.Print();
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
    PlasticityFunction(epsbar, c, H);
    
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
    sigma_projected = sigma_trial;
    NRvector<T>  eigenvalues = sigma_projected;
    const Doub sinphi = sin(fPhi);
    const Doub sinpsi = sin(fPsi);
    const Doub cosphi = cos(fPhi);
    const Doub sinphi2 = sinphi*sinphi;
    const Doub cosphi2 = 1. - sinphi2;
    NRvector<T> gamma(2, 0.), phi(2, 0.), sigma_bar(2, 0.), ab(2, 0.);
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    
    NRvector<T> phival(2, 0.);
    NRvector<NRvector<T>> jac(2), jac_inv(2);
    for (int i = 0; i < 2; i++) {
        jac[i].assign(2, 0.);
        jac_inv[i].assign(2, 0.);
    }
    
    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * T(sinphi);
    sigma_bar[1] = eigenvalues[1] - eigenvalues[2]+(eigenvalues[1] + eigenvalues[2]) * T(sinphi);
    T c, H;
    ab[0] = T(4. * fmu*(1 + sinphi * sinpsi / 3.) + 4. * fK* sinphi * sinpsi);
    ab[1] = T(2. * fmu*(1. - sinphi - sinpsi - sinphi * sinpsi / 3.) + 4. * fK * sinphi * sinpsi);
    T epsbar = T(fEpsPlasticBar) + (gamma[0] + gamma[1]) * T(2. * cosphi);
    PlasticityFunction(epsbar, c, H);
    
    phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
    phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;
    T res = (fabs(phival[0]) + fabs(phival[1]));
    int n_iterations = 30;  
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {
        
        jac[0][0] = -ab[0] - T(4. * cosphi2) * H;
        jac[1][0] = -ab[1] - T(4. * cosphi2) * H;
        jac[0][1] = -ab[1] - T(4. * cosphi2) * H;
        jac[1][1] = -ab[0] - T(4. * cosphi2) * H;
        
        T det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        

        if(fabs(det_jac)<1.e-8){
            std::cerr << "MohrCoulomb:: Singular jacobian." << std::endl;
            DebugStop();
        }

        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;
        
        gamma[0] -= (jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1]);
        gamma[1] -= (jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1]);
        
        epsbar = T(fEpsPlasticBar)+(gamma[0] + gamma[1]) * T(2. * cosphi);
        PlasticityFunction(epsbar, c, H);
        
        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;
        res = (fabs(phi[0]) + fabs(phi[1]));
        
        if (fabs(res)<1.e-8) {
            break;
        }
    }
    
    if (i == n_iterations) {
        DebugStop();
    }

    memory.fGamma[0] = gamma[0];
    memory.fGamma[1] = gamma[1];
    eigenvalues[0] += -T(2. * fmu*(1 + sinpsi / 3.) + 2. * fK * sinpsi) * gamma[0] + T((4. * fmu / 3. - 2. * fK) * sinpsi) * gamma[1];
    eigenvalues[1] += T((4. * fmu / 3. - fK*2.) * sinpsi) * gamma[0] - T(2. * fmu*(1. + sinpsi / 3.) + 2. * fK * sinpsi) * gamma[1];
    eigenvalues[2] += T(2. * fmu*(1 - sinpsi / 3.) - 2. * fK * sinpsi)*(gamma[0] + gamma[1]);
    sigma_projected = eigenvalues;
    epsbarnew = epsbar;

    bool check_validity_Q = ( eigenvalues[0]> eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1])<1.e-8) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1]-eigenvalues[2])<1.e-8);
    return (check_validity_Q);
}

void mohrcoulomb::ComputeLeftEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    const Doub sin_phi = sin(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub G = fmu, K = fK;
    const Doub a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    const Doub b = 2.0 * G * (1.0 - sin_phi - sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    
    Doub epsbar = epsbarp;
    Doub c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.assign(3, 3,0.);
    
    // First column
    tang[0][0] = (-3*b*b + 3*a*(a - 2*G*(1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    tang[1][0] = (2*(1 + sin_phi)*(a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi + b*G*(3 + sin_psi)))/
    (3.*(a*a - b*b));
    tang[2][0] = (-2*(1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));
    
    // Second column
    tang[0][1] = (2*(1 + sin_phi)*(a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi + b*G*(3 + sin_psi)))/
    (3.*(a*a - b*b));
    tang[1][1] = (-3*b*b + 3*a*(a - 2*G*(1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    tang[2][1] = (-2*(1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));
    
    // Third column
    tang[0][2] = (2*(-1 + sin_phi)*(G*(-3 + sin_psi) - 6*K*sin_psi))/(3.*(a + b));
    tang[1][2] = (2*(-1 + sin_phi)*(G*(-3 + sin_psi) - 6*K*sin_psi))/(3.*(a + b));
    tang[2][2] = (3*a + 3*b - 4*(-1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));
    
}

template<class T>
bool mohrcoulomb::ReturnMapRightEdge(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const
{
    sigma_projected = sigma_trial;
    NRvector<T> eigenvalues = sigma_projected;
    const Doub sinphi = sin(fPhi);
    const Doub sinpsi = sin(fPsi);
    const Doub cosphi = cos(fPhi);
    const Doub sinphi2 = sinphi*sinphi;
    const Doub cosphi2 = 1. - sinphi2;
    const Doub KV = fK;
    const Doub GV = fmu;
    
    NRvector<T> gamma(2, 0.), phi(2, 0.), sigma_bar(2, 0.), ab(2, 0.);
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    NRvector<T> phival(2, 0.);
    NRvector<NRvector<T>> jac(2), jac_inv(2);
    for (int i = 0; i < 2; i++) {
        jac[i].assign(2, 0.);
        jac_inv[i].assign(2, 0.);
    }
    
    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * T(sinphi);
    sigma_bar[1] = eigenvalues[0] - eigenvalues[1]+(eigenvalues[0] + eigenvalues[1]) * T(sinphi);
    T c, H;
    ab[0] = T(4. * GV * (1 + sinphi * sinpsi / 3.) + 4. * KV * sinphi * sinpsi);
    ab[1] = T(2. * GV * (1. + sinphi + sinpsi - sinphi * sinpsi / 3.) + 4. * KV * sinphi * sinpsi);
    T epsbar = T(fEpsPlasticBar)+(gamma[0] + gamma[1]) * T(2. * cosphi);
    PlasticityFunction(epsbar, c, H);
    
    phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
    phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;

    T res = (fabs(phival[0]) + fabs(phival[1]));
    int n_iterations = 30; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {

        jac[0][0] = -ab[0] - T(4. * cosphi2) * H;
        jac[1][0] = -ab[1] - T(4. * cosphi2) * H;
        jac[0][1] = -ab[1] - T(4. * cosphi2) * H;
        jac[1][1] = -ab[0] - T(4. * cosphi2) * H;
        
        T det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

        if(fabs(det_jac)<1.e-8){
            std::cerr << "MohrCoulombPV:: Singular jacobian." << std::endl;
            DebugStop();
        }

        
        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;
        
        gamma[0] -= (jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1]);
        gamma[1] -= (jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1]);
        
        epsbar = T(fEpsPlasticBar)+(gamma[0] + gamma[1]) * T(2. * cosphi);
        PlasticityFunction(epsbar, c, H);

        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T(2. * cosphi) * c;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T(2. * cosphi) * c;
        phival[0] = phi[0];
        phival[1] = phi[1];
        res = (fabs(phival[0]) + fabs(phival[1]));
        
        
        if (fabs(res)<1.e-8) {
            break;
        }
    }
    
    if (i == n_iterations) {
        DebugStop();
    }


    memory.fGamma[0] = gamma[0];
    memory.fGamma[1] = gamma[1];

    eigenvalues[0] -= T(2. * GV * (1 + sinpsi / 3.) + 2. * KV * sinpsi)*(gamma[0] + gamma[1]);
    eigenvalues[1] += T((4. * GV / 3. - KV * 2.) * sinpsi) * gamma[0] + T(2. * GV * (1. - sinpsi / 3.) - 2. * KV * sinpsi) * gamma[1];
    eigenvalues[2] += T(2. * GV * (1 - sinpsi / 3.) - 2. * KV * sinpsi) * gamma[0] + T((4. * GV / 3. - 2. * KV) * sinpsi) * gamma[1];
    sigma_projected = eigenvalues;
    epsbarnew = epsbar;

    bool check_validity_Q = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1])<1.e-8) && ( eigenvalues[1]  >  eigenvalues[2]  || fabs(eigenvalues[1]-eigenvalues[2])<1.e-8);
    return (check_validity_Q);
}

void mohrcoulomb::ComputeRightEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    const Doub sin_phi = sin(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub G = fmu, K = fK;
    const Doub a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    const Doub b = 2.0 * G * (1.0 + sin_phi + sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;

    Doub epsbar = epsbarp;
    Doub c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.assign(3, 3,0.);
    
    // First column
    tang[0][0] = (3.0*a + 3.0*b - 4.0*(1.0 + sin_phi)*(3.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
    tang[1][0] = (2.0*(1.0 + sin_phi)*(-6.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
    tang[2][0] = (2.0*(1.0 + sin_phi)*(-6.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
    
    // Second column
    tang[0][1] = (-2*(-1 + sin_phi)*(3*K*sin_psi + G*(3 + sin_psi)))/(3.*(a + b));
    tang[1][1] = (-3*b*b + 3*a*(a + 2*G*(-1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
    tang[2][1] = (2*(-1 + sin_phi)*(b*G*(-3 + sin_psi) + a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi))/
    (3.*(a*a - b*b));
    
    // Third column
    tang[0][2] = (-2*(-1 + sin_phi)*(3*K*sin_psi + G*(3 + sin_psi)))/(3.*(a + b));
    tang[1][2] = (2*(-1 + sin_phi)*(b*G*(-3 + sin_psi) + a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi))/
    (3.*(a*a - b*b));
    tang[2][2] = (-3*b*b + 3*a*(a + 2*G*(-1 + sin_phi)) -
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
}

template<class T>
bool mohrcoulomb::ReturnMapApex(const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
            TComputeSequence &memory, Doub &epsbarnew) const
{
    const Doub K = fK;
    const Doub sinpsi = sin(fPsi);
    const Doub cosphi = cos(fPhi);
    const Doub cotphi = 1. / tan(fPhi);
    T ptrnp1 = 0.;
    for (int i = 0; i < 3; i++) {
        ptrnp1 += T(sigma_trial[i]);
    }
    ptrnp1 /= 3.;
    T DEpsPV = 0.;
    T epsbarnp1 = T(fEpsPlasticBar);
    T c, H;
    PlasticityFunction(epsbarnp1, c, H);

    T alpha = cos(fPhi) / sin(fPsi);
    Doub tol = 1.e-8;

    T res = c * cotphi - ptrnp1;
    T pnp1;
    
    int n_iterations = 30; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for (i = 0; i < n_iterations; i++) {
        const T jac = H * T(cosphi * cotphi) / T(sinpsi) + T(K);
        DEpsPV -= res / jac;

        epsbarnp1 = T(fEpsPlasticBar) + T(alpha) * DEpsPV;
        pnp1 = ptrnp1 - T(K) * DEpsPV;
        PlasticityFunction(epsbarnp1, c, H);
        res = c * cotphi - pnp1;
        
        if (fabs(res)<1.e-8) {
            break;
        }
    }
    
    if (i == n_iterations) {
        DebugStop();
    }

    epsbarnew = epsbarnp1;
    for (int i = 0; i < 3; i++) {
        sigma_projected[i] = pnp1;
    }
    return true;
}
    

void mohrcoulomb::ComputeApexGradient(NRmatrix<Doub> & gradient, Doub & eps_bar_p) const
{
    Doub c, H;
    const Doub cosphi = cos(fPhi);
    const Doub sinpsi = sin(fPsi);
    const Doub cotphi = 1. / tan(fPhi);
    const Doub K = fK;
    const Doub alpha = cosphi / sinpsi;
    PlasticityFunction(eps_bar_p, c, H);
    const Doub num = H * alpha * cotphi / K;
    const Doub denom = 1. + num;
    const Doub dpdptr = num / denom;
    const Doub dsigdsigtr = dpdptr / 3.;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            gradient[i][j] = dsigdsigtr;
        }
    }
}
    

void mohrcoulomb::ProjectSigma(const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub>  gradient)
{

    
    TComputeSequence memory;
    this->SetEpsBar(k_prev);
    Doub epsbartemp = k_prev; // it will be defined by the correct returnmap
    
    bool check_validity_Q;

    // Check if we are in the correct sextant
    check_validity_Q = (sigma_trial[0] > sigma_trial[1] || fabs(sigma_trial[0]-sigma_trial[1])<1.e-8) && (sigma_trial[1] > sigma_trial[2] || fabs(sigma_trial[1]-sigma_trial[2])<1.e-8);
    if (!check_validity_Q) {
        DebugStop();
    }

//22652.243938097836
    Doub phi = PhiPlane<Doub>(sigma_trial);
    bool elastic_update_Q = fabs(phi)<1.e-8 || phi < 0.0;
    if (elastic_update_Q) {
        m_type = 0; // Elastic behavior
        memory.fWhichPlane = TComputeSequence::EElastic;
        memory.fGamma.resize(0);
        sigma = sigma_trial;
        
        
            for(Int i =0;i<3;i++)
            {
                gradient[i][i] = 1.;
            }
        
        return;
    }
    
    m_type = 1; // failure behavior
    NRvector<Doub>sigma_projected;
    memory.fGamma.resize(1);
    memory.fGamma[0] = 0.;
    check_validity_Q = this->ReturnMapPlane<Doub>(sigma_trial, sigma_projected, memory, epsbartemp);
    if (check_validity_Q) {
        k_proj = epsbartemp;
        this->SetEpsBar(k_proj);
        sigma = sigma_projected;
        memory.fWhichPlane = TComputeSequence::EMainPlane;
        
        ComputePlaneTangent(gradient, epsbartemp);
        
    } else {
        memory.fGamma.resize(2);
        memory.fGamma[0] = 0.;
        memory.fGamma[1] = 0.;
        bool IsEdge = false;

        const Doub sinpsi = sin(fPsi);
        Doub val = (1 - sinpsi) * sigma_trial[0] - 2. * sigma_trial[1] + (1 + sinpsi) * sigma_trial[2];
        if (val > 0.) {
            IsEdge = this->ReturnMapRightEdge<Doub>(sigma_trial, sigma_projected, memory, epsbartemp);
            memory.fWhichPlane = TComputeSequence::ERightEdge;
            
            ComputeRightEdgeTangent(gradient, epsbartemp);
            
        } else {
            IsEdge = this->ReturnMapLeftEdge<Doub>(sigma_trial, sigma_projected, memory, epsbartemp);
            memory.fWhichPlane = TComputeSequence::ELeftEdge;
            
            ComputeLeftEdgeTangent(gradient, epsbartemp);

        }
        if (!IsEdge) {
            m_type = -1; // Tensile behavior
            this->ReturnMapApex(sigma_trial, sigma_projected, memory, epsbartemp);
            memory.fWhichPlane = TComputeSequence::EApex;
            
            ComputeApexGradient(gradient, epsbartemp);
        }

        k_proj = epsbartemp;
        this->SetEpsBar(k_proj);
        sigma = sigma_projected;
    }
}

/*
void mohrcoulomb::Phi(NRvector<Doub> sig_vec, Doub alpha, NRvector<Doub> &phi)const
{
    phi.resize(3);
    for (int i = 0; i < 3; i++) phi[i] = 0;
    phi[0] = PhiPlane(sig_vec);
    phi[2] = PhiPlane(sig_vec); // Consistency with two surfaces models
}
*/
