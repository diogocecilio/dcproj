#include "mohrcoulomb.h"

mohrcoulomb::mohrcoulomb()
{
    //DebugStop();
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
    
    NRvector<Doub>  phit(3);
    Doub alpha = 0.0;
    Doub tol = 1.0e-10;
    Phi(stress_p, alpha, phit);
    bool Is_valid_stress_on_cap_Q =  fabs(phit[0]) < tol || phit[0] < 0.0;
    
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
        if (fabs(phi)<ftol) {
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
    //1sigma_projected.Print();
    epsbarnew = epsbar;
    
    bool check_validity_Q = (eigenvalues[0] > eigenvalues[1]|| fabs(eigenvalues[0]-eigenvalues[1])<ftol && (eigenvalues[1] > eigenvalues[2]) || (eigenvalues[1]-eigenvalues[2])<ftol);
    return (check_validity_Q);   
}

void mohrcoulomb::ComputePlaneTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    const Doub sin_phi = sin(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub cos_phi = cos(fPhi);
    const Doub G = fmu, K = fK;
    const Doub denominator = 6.0 * G + 2.0 * (G + 3.0 * K) * sin_phi * sin_psi;
    
    Doub epsbar = epsbarp;
    Doub c=fc, H=0;
    PlasticityFunction(epsbar, c, H);
    
    tang.assign(3, 3,0.);
    
    // First column
   /* tang[0][0] = (sin_phi - 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    tang[1][0] = (2.0 * G - 3.0 * K) * (sin_phi + 1.0) * sin_psi / denominator;
    tang[2][0] = -(sin_phi + 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    
    // Second column
    tang[0][1] = 0.0;
    tang[1][1] = 1.0;
    tang[2][1] = 0.0;
    
    // Third column
    tang[0][2] = -(sin_phi - 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
    tang[1][2] = (2.0 * G - 3.0 * K) * (sin_phi - 1.0) * sin_psi / denominator;
    tang[2][2] = (sin_phi + 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;*/
    
    Doub R4G =  4.*fmu;
    Doub R2G = 2. *fmu;
    Doub R1 =1.;
    Doub R2 =2.;
    Doub R3=3.;
    Doub SINPHI = sin_phi;
    Doub SINPSI = sin_psi;
    Doub R1D3 = 1./3.;
    Doub SPHSPS = sin_phi*sin_psi;
    Doub R4=4.;
    Doub BULK =fK;
    Doub CONSTB=R2G*(R1-SINPHI-SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS;
    Doub  R2CPHI= 2.*cos_phi;
    Doub R4C2PH = R2CPHI*R2CPHI;
    Doub R2GD3=R2G*R1D3;
    Doub R4GD3=R4G*R1D3;
    
    Doub CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS;
    Doub R2BULK = 2*fK;
    Doub FACTA=R4C2PH*H;
    Doub DRVAA=-CONSTA-FACTA;
    Doub DRVAB=-CONSTB-FACTA;
    Doub DRVBA=-CONSTB-FACTA;
    Doub DRVBB=-CONSTA-FACTA;
    Doub AUX1=R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI;
    Doub AUX2=(R4GD3-R2BULK)*SINPSI;
    Doub AUX3=R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI;
    Doub R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA);
    Int II=0;
    Int JJ=1;
    Int MM=2;
    Doub R2D3=2./3.;
    
 
    CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS;
    
    Doub DENOM=-CONSTA-R4C2PH*H;
    Doub B1=(R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI)/DENOM;
    Doub B2=(R4G*R1D3-R2BULK)*SINPSI/DENOM;
    Doub B3=(R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI)/DENOM;
    
          
    tang[II][II]=R2G*(R2D3+B1*(R1+R1D3*SINPHI))+BULK*(R1+R2*B1*SINPHI);
    
    tang[II][MM]=R1D3*(R3*BULK-R2G)*(R1+R2*B1*SINPHI);
          
    tang[II][JJ]=R2G*(-R1D3-B1*(R1-R1D3*SINPHI))+BULK*(R1+R2*B1*SINPHI);
    
    tang[MM][II]=R2G*(-R1D3-B2*(R1+R1D3*SINPHI))+BULK*(R1-R2*B2*SINPHI);
          
    tang[MM][MM]=R4G*R1D3*(R1+B2*SINPHI)+BULK*(R1-R2*B2*SINPHI);
    
    tang[MM][JJ]=R2G*(-R1D3+B2*(R1-R1D3*SINPHI))+BULK*(R1-R2*B2*SINPHI);
    
    tang[JJ][II]=R2G*(-R1D3-B3*(R1+R1D3*SINPHI))+BULK*(R1-R2*B3*SINPHI);
    
    tang[JJ][MM]=R1D3*(R3*BULK-R2G)*(R1-R2*B3*SINPHI);
    
    tang[JJ][JJ]=R2G*(R2D3+B3*(R1-R1D3*SINPHI))+BULK*(R1-R2*B3*SINPHI);
    
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
        

        if(fabs(det_jac)<ftol){
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
        
        if (fabs(res)<ftol) {
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

    bool check_validity_Q = ( eigenvalues[0]> eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1])<ftol) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1]-eigenvalues[2])<ftol);
    return (check_validity_Q);
}

void mohrcoulomb::ComputeLeftEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    const Doub sin_phi = sin(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub cos_phi = cos(fPhi);
    const Doub G = fmu, K = fK;
    const Doub a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    const Doub b = 2.0 * G * (1.0 - sin_phi - sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    
    Doub epsbar = epsbarp;
    Doub c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.assign(3, 3,0.);
  /* 
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
    */
  
    Doub R4G =  4.*fmu;
    Doub R2G = 2. *fmu;
    Doub R1 =1.;
    Doub SINPHI = sin_phi;
    Doub SINPSI = sin_psi;
    Doub R1D3 = 1./3.;
    Doub SPHSPS = sin_phi*sin_psi;
    Doub R4=4.;
    Doub BULK =fK;
    Doub CONSTB=R2G*(R1-SINPHI-SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS;
    Doub  R2CPHI= 2.*cos_phi;
    Doub R4C2PH = R2CPHI*R2CPHI;
    Doub R2GD3=R2G*R1D3;
    Doub R4GD3=R4G*R1D3;
    
    Doub CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS;
    Doub R2BULK = 2*fK;
    Doub FACTA=R4C2PH*H;
    Doub DRVAA=-CONSTA-FACTA;
    Doub DRVAB=-CONSTB-FACTA;
    Doub DRVBA=-CONSTB-FACTA;
    Doub DRVBB=-CONSTA-FACTA;
    Doub AUX1=R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI;
    Doub AUX2=(R4GD3-R2BULK)*SINPSI;
    Doub AUX3=R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI;
    Doub R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA);
    Int II=0;
    Int JJ=1;
    Int MM=2;
    
    tang[II][II] = BULK+R4GD3+(AUX1*(((R2BULK*(DRVBB-DRVAB)+(DRVAB*R4GD3+DRVBB*R2GD3))*SINPHI)+DRVBB*R2G)+AUX2*(((R2BULK*(DRVBA-DRVAA)+(DRVAA*R4GD3+DRVBA*R2GD3))*SINPHI)+DRVBA*R2G))*R1DDET;
    
    tang[II][MM]=BULK-R2GD3+(AUX1*(((R2BULK*(DRVBB-DRVAB)-(DRVAB*R2GD3+DRVBB*R4GD3))*SINPHI)-DRVAB*R2G)+AUX2*(((R2BULK*(DRVBA-DRVAA)-(DRVAA*R2GD3+DRVBA*R4GD3))*SINPHI)-DRVAA*R2G))*R1DDET;
                         
    tang[II][JJ]=BULK-R2GD3+((AUX1*(DRVBB-DRVAB)+AUX2*(DRVBA-DRVAA))*(((R2BULK+R2GD3)*SINPHI)-R2G))*R1DDET;
    
    tang[MM][II]=BULK-R2GD3+(AUX1*(((R2BULK*(DRVAA-DRVBA)-(DRVAA*R4GD3+DRVBA*R2GD3))*SINPHI)-DRVBA*R2G)+AUX2*(((R2BULK*(DRVAB-DRVBB)-(DRVAB*R4GD3+DRVBB*R2GD3))*SINPHI)-DRVBB*R2G))*R1DDET;
    
    tang[MM][MM]=BULK+R4GD3+(AUX1*(((R2BULK*(DRVAA-DRVBA)+(DRVAA*R2GD3+DRVBA*R4GD3))*SINPHI)+DRVAA*R2G)+AUX2*(((R2BULK*(DRVAB-DRVBB)+(DRVAB*R2GD3+DRVBB*R4GD3))*SINPHI)+DRVAB*R2G))*R1DDET;
    
    tang[MM][JJ] = BULK-R2GD3+((AUX1*(DRVAA-DRVBA)+AUX2*(DRVAB-DRVBB))*(((R2BULK+R2GD3)*SINPHI)-R2G))*R1DDET;
    
    tang[JJ][II]=BULK-R2GD3+(AUX3*(((R2BULK*(DRVAB-DRVBB-DRVAA+DRVBA)+(DRVAA-DRVAB)*R4GD3+(DRVBA-DRVBB)*R2GD3)*SINPHI)+(DRVBA-DRVBB)*R2G))*R1DDET;
    
    tang[JJ][MM]=BULK-R2GD3+(AUX3*(((R2BULK*(DRVAB-DRVBB-DRVAA+DRVBA)+(DRVAB-DRVAA)*R2GD3+(DRVBB-DRVBA)*R4GD3)*SINPHI)+(DRVAB-DRVAA)*R2G))*R1DDET;
            
    tang[JJ][JJ]=BULK+R4GD3+(AUX3*(DRVAB-DRVBB-DRVAA+DRVBA)*(((R2BULK+R2GD3)*SINPHI)-R2G))*R1DDET;
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

        if(fabs(det_jac)<ftol){
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
        
        
        if (fabs(res)<ftol) {
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

    bool check_validity_Q = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1])<ftol) && ( eigenvalues[1]  >  eigenvalues[2]  || fabs(eigenvalues[1]-eigenvalues[2])<ftol);
    return (check_validity_Q);
}

void mohrcoulomb::ComputeRightEdgeTangent(NRmatrix<Doub> &tang, Doub &epsbarp) const
{
    const Doub sin_phi = sin(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub cos_phi = cos(fPhi);
    const Doub G = fmu, K = fK;
    const Doub a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
    const Doub b = 2.0 * G * (1.0 + sin_phi + sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;

    Doub epsbar = epsbarp;
    Doub c, H;
    PlasticityFunction(epsbar, c, H);
    
    tang.assign(3, 3,0.);
    
    // First column
   /* tang[0][0] = (3.0*a + 3.0*b - 4.0*(1.0 + sin_phi)*(3.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
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
                  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));*/
                  
                  
    Doub R4G =  4.*fmu;
    Doub R2G = 2. *fmu;
    Doub R1 =1.;
    Doub SINPHI = sin_phi;
    Doub SINPSI = sin_psi;
    Doub R1D3 = 1./3.;
    Doub SPHSPS = sin_phi*sin_psi;
    Doub R4=4.;
    Doub BULK =fK;
    Doub CONSTB=R2G*(R1+SINPHI+SINPSI-R1D3*SPHSPS)+R4*BULK*SPHSPS;
    Doub  R2CPHI= 2.*cos_phi;
    Doub R4C2PH = R2CPHI*R2CPHI;
    Doub R2GD3=R2G*R1D3;
    Doub R4GD3=R4G*R1D3;
    
    Doub CONSTA=R4G*(R1+R1D3*SPHSPS)+R4*BULK*SPHSPS;
    Doub R2BULK = 2*fK;
    Doub FACTA=R4C2PH*H;
    Doub DRVAA=-CONSTA-FACTA;
    Doub DRVAB=-CONSTB-FACTA;
    Doub DRVBA=-CONSTB-FACTA;
    Doub DRVBB=-CONSTA-FACTA;
    Doub AUX1=R2G*(R1+R1D3*SINPSI)+R2BULK*SINPSI;
    Doub AUX2=(R4GD3-R2BULK)*SINPSI;
    Doub AUX3=R2G*(R1-R1D3*SINPSI)-R2BULK*SINPSI;
    Doub R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA);
    Int II=0;
    Int JJ=1;
    Int MM=2;
    
    tang[II][II]=BULK+R4GD3+AUX1*(-DRVAB+DRVBB+DRVAA-DRVBA)*
                         (R2G+(R2BULK+R2GD3)*SINPHI)*R1DDET;
                         
    tang[II][MM]=BULK-R2GD3+AUX1*(R2G*(DRVAB-DRVAA)+((-DRVAB+DRVBB+DRVAA-DRVBA)*(R2BULK+R2GD3)+(DRVBA-DRVBB)*R2G)*SINPHI)*R1DDET;
                         
    tang[II][JJ]=BULK-R2GD3+AUX1*(R2G*(DRVBA-DRVBB)+((-DRVAB+DRVBB+DRVAA-DRVBA)*(R2BULK+R2GD3)+(DRVAB-DRVAA)*R2G)*SINPHI)*R1DDET;
                         
    tang[MM][II]=BULK-R2GD3+(AUX2*(DRVAB-DRVBB)+AUX3*(DRVBA-DRVAA))*(R2G+(R2BULK+R2GD3)*SINPHI)*R1DDET;
    
    tang[MM][MM]=BULK+R4GD3+(AUX2*((R2BULK*(DRVAB-DRVBB)+(DRVAB*R2GD3+DRVBB*R4GD3))*SINPHI-DRVAB*R2G)+AUX3*(DRVAA*R2G+(R2BULK*(DRVBA-DRVAA)-
        (DRVAA*R2GD3+DRVBA*R4GD3))*SINPHI))*R1DDET;
        
    tang[MM][JJ]=BULK-R2GD3+(AUX2*((R2BULK*(DRVAB-DRVBB)-(DRVBB*R2GD3+DRVAB*R4GD3))*SINPHI+DRVBB*R2G)+AUX3*((R2BULK*(DRVBA-DRVAA)+(DRVAA*R4GD3+DRVBA*R2GD3))*SINPHI-DRVBA*R2G))*R1DDET;
                         
    tang[JJ][II]=BULK-R2GD3+((AUX2*(DRVBA-DRVAA)+AUX3*(DRVAB-DRVBB))*((R2BULK+R2GD3)*SINPHI+R2G))*R1DDET;
    
    tang[JJ][MM]=BULK-R2GD3+(AUX2*(((R2BULK*(DRVBA-DRVAA)-(DRVBA*R4GD3+DRVAA*R2GD3))*SINPHI)+DRVAA*R2G)+AUX3*(((R2BULK*(DRVAB-DRVBB)+(DRVAB*R2GD3+DRVBB*R4GD3))*SINPHI)-DRVAB*R2G))*R1DDET;
    
    tang[JJ][JJ]=BULK+R4GD3+(AUX2*(((R2BULK*(DRVBA-DRVAA)+(DRVAA*R4GD3+DRVBA*R2GD3))*SINPHI)-DRVBA*R2G)+AUX3*(((R2BULK*(DRVAB-DRVBB)-(DRVAB*R4GD3+DRVBB*R2GD3))*SINPHI)+DRVBB*R2G))*R1DDET;
                  
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
    Doub tol = ftol;

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
        
        if (fabs(res)<ftol) {
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
    const Doub cos_phi = cos(fPhi);
    const Doub sin_psi = sin(fPsi);
    const Doub sin_phi = sin(fPhi);
    const Doub cotphi = 1. / tan(fPhi);
    const Doub K = fK;
    const Doub alpha = cos_phi / sin_psi;
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
    
    Doub R4G =  4.*fmu;
    Doub R2G = 2. *fmu;
    Doub R1 =1.;

    Doub BULK =fK;
    Doub COSPHI = cos_phi;
    Doub SINPHI = sin_phi;
    Doub SINPSI = sin_psi;
    
    
        Doub COTPHI=COSPHI/SINPHI;
        Doub DSIDEJ=BULK*(R1-(BULK/(BULK+H*COTPHI*COSPHI/SINPSI)));
          
        Int II =0,JJ=1,MM=2;
          
        gradient[II][II]=DSIDEJ;
        gradient[II][MM]=DSIDEJ;
        gradient[II][JJ]=DSIDEJ;
        gradient[MM][II]=DSIDEJ;
        gradient[MM][MM]=DSIDEJ;
        gradient[MM][JJ]=DSIDEJ;
        gradient[JJ][II]=DSIDEJ;
        gradient[JJ][MM]=DSIDEJ;
        gradient[JJ][JJ]=DSIDEJ;
}
    

void mohrcoulomb::ProjectSigma(const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub>  &gradient)
{

    gradient.assign(3,3,0.);
    TComputeSequence memory;
    this->SetEpsBar(k_prev);
    Doub epsbartemp = k_prev; // it will be defined by the correct returnmap
    
    bool check_validity_Q;

    // Check if we are in the correct sextant
    check_validity_Q = (sigma_trial[0] > sigma_trial[1] || fabs(sigma_trial[0]-sigma_trial[1])<ftol) && (sigma_trial[1] > sigma_trial[2] || fabs(sigma_trial[1]-sigma_trial[2])<ftol);
    if (!check_validity_Q) {
        DebugStop();
    }

//22652.243938097836
    Doub phi = PhiPlane<Doub>(sigma_trial);
    bool elastic_update_Q = fabs(phi)<ftol || phi < 0.0;
    if (elastic_update_Q) {
        m_type = 0; // Elastic behavior
        memory.fWhichPlane = TComputeSequence::EElastic;
        memory.fGamma.resize(0);
        sigma = sigma_trial;
        
        NRmatrix<Doub> C = GetElasticMatrix();
        gradient[0][0] = C[0][0];gradient[0][1] = C[0][1];gradient[0][2] = C[0][5];
        gradient[1][0] = C[1][0];gradient[1][1] = C[1][1];gradient[1][2] = C[1][5];
        gradient[2][0] = C[5][0];gradient[2][1] = C[5][1];gradient[2][2] = C[5][5];
        
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

void mohrcoulomb::updateatributes(NRvector<MatDoub> mult)
{
	Doub newcoesion = 0.;
	Doub newphi = 0.;
	//Doub multcoes = 0., multphi = 0.;
	fyoung0 = fyoung;
	fnu0 = fnu;
	fc0 = fc;
	fPhi0 = fPhi;
    fPsi0 = fPsi;
	//Doub newyoung = fyoung + mult*fyoung;
	//Doub newcoesion = fcoesion + mult[0][0][0]*fcoesion;
	//Doub newphi = fphi + mult[1][0][0] * fphi;

	newcoesion = mult[0][0][0];
	newphi = mult[1][0][0];
	SetUp(fPhi, fPsi, fc, fyoung,fnu);
}

void mohrcoulomb::closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma)
{
	NRtensor<Doub>  epse = epst - epsp;
	NRmatrix<Doub>  C = GetElasticMatrix();

	NRmatrix<Doub>  tempepsemat, stresstrial, Dept;
	epse.FromTensorToNRmatrix(tempepsemat);

	C.Mult(tempepsemat, stresstrial);

	NRtensor<Doub>  stresstrialtensor;
	epse.FromNRmatrixToTensor(stresstrial, stresstrialtensor);

    NRmatrix<Doub>  pt, vec;
    NRvector<Doub> sigma_trial(3,0.),sigma,gradient(3,3);
    
    stresstrialtensor.EigenSystem(pt, vec);
    
    sigma_trial[0] = pt[0][0];
    sigma_trial[1] =  pt[0][1];
    sigma_trial[2] = pt[0][2];
			
    

    Doub k_prev=0.;
    Doub k_proj=0.;
    Int  m_type;
    ProjectSigma( sigma_trial, k_prev, sigma, k_proj,m_type, Dep);

}
