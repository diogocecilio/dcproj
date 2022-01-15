#include "mohrcoulomb.h"

mohrcoulomb::mohrcoulomb()
{
    //DebugStop();
}
mohrcoulomb::mohrcoulomb ( Doub young, Doub nu, Doub c, Doub Phi,Doub Psi )
{
    fPhi = Phi;
    fPsi = Psi;
    fc = c;
    fyoung = young;
    fnu = nu;
    SetUp ( young,  nu,c,Phi,Psi );
}
mohrcoulomb::mohrcoulomb ( const mohrcoulomb &cp )
{
    DebugStop();
}

Doub mohrcoulomb::InitialDamage ( const NRvector<Doub> &stress_p ) const
{

    std::cout << "TPZYCMohrCoulombPV::There is no damage variable for this model at the current time." << std::endl;

    NRvector<Doub>  phit ( 3 );
    Doub alpha = 0.0;
    Doub tol = 1.0e-10;
    Phi ( stress_p, alpha, phit );
    bool Is_valid_stress_on_cap_Q =  fabs ( phit[0] ) < tol || phit[0] < 0.0;

    if ( !Is_valid_stress_on_cap_Q ) {
        std::cerr << "MohrCoulombPV::Invalid stress state." << std::endl;
        DebugStop();
    }

    return 0.0;
}

void mohrcoulomb::Phi ( NRvector<Doub>sig_vec, Doub alpha, NRvector<Doub> &phi ) const
{
    phi.assign ( 3,0. );
    phi[0] = PhiPlane ( sig_vec );
}

template <class T>
void mohrcoulomb::PlasticityFunction ( const T epsp, T &c, T &H ) const
{
    c = fc; // c(epsp)
    H = 0.00001; // dc(epsp)/depsp
}

template<class T>
NRvector<T> mohrcoulomb::SigmaElastPV ( const NRvector<T> &deform ) const
{
    T trace = deform[0] + deform[1] + deform[2];
    NRvector<T> sigma ( 3, 0. );


    sigma = trace * flambda + 2 * fmu * deform[0];
    sigma = trace * flambda+ 2 *fmu * deform[1];
    sigma = trace * flambda + 2 * fmu * deform[2];

    return sigma;
}

template<class T>
T mohrcoulomb::PhiPlane ( const NRvector<T> &sigma ) const
{
    const Doub sinphi = sin ( fPhi );
    const Doub cosphi = cos ( fPhi );
    Doub temp =sigma[0] - sigma[2] + ( sigma[0] + sigma[2] ) * sinphi;
    return  temp - 2. * fc*cosphi;
}


template<class T>
bool mohrcoulomb::ReturnMapPlane ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                                   TComputeSequence &memory, Doub &epsbarnew ) const
{
    sigma_projected = sigma_trial;
    NRvector<T> eigenvalues = sigma_projected;
    const Doub sinphi = sin ( fPhi );
    const Doub sinpsi = sin ( fPsi );
    const Doub cosphi = cos ( fPhi );
    const Doub sinphi2 = sinphi*sinphi;
    const Doub cosphi2 = 1. - sinphi2;
    const Doub constA = 4. * fmu * ( 1. + sinphi * sinpsi / 3. ) + 4. * fK * sinphi*sinpsi;
    T c, H;
    T epsbar = T ( fEpsPlasticBar + memory.fGamma[0]*2. * cosphi );
    PlasticityFunction ( epsbar, c, H );

    T phi = eigenvalues[0] - eigenvalues[2]+ ( eigenvalues[0] + eigenvalues[2] ) * sinphi - 2. * fc*cosphi;
    T gamma = memory.fGamma[0];
    int n_iterations = 300;

    int i;
    bool stop_criterion;
    for ( i = 0; i < n_iterations; i++ ) {
        T jac = -constA - T ( 4. * cosphi2 ) * H;
        T delta_gamma = - phi / jac;
        gamma += delta_gamma;
        phi = eigenvalues[0] - eigenvalues[2]+ ( eigenvalues[0] + eigenvalues[2] ) * sinphi - 2. * fc * cosphi - constA * gamma;
        if ( fabs ( phi ) <ftol ) {
            break;
        }
    }

    if ( i == n_iterations ) {
        //DebugStop();
    }

    memory.fGamma[0] = gamma;
    eigenvalues[0] -= T ( 2. * fmu* ( 1. + sinpsi / 3. ) + 2. * fK * sinpsi ) * gamma;
    eigenvalues[1] += T ( ( 4. * fmu / 3. - fK*2. ) * sinpsi ) * gamma;
    eigenvalues[2] += T ( 2. * fmu* ( 1 - sinpsi / 3. ) - 2. * fK * sinpsi ) * gamma;
    sigma_projected = eigenvalues;
    //1sigma_projected.Print();
    epsbarnew = epsbar;

    bool check_validity_Q = ( eigenvalues[0] > eigenvalues[1] || fabs ( eigenvalues[0]-eigenvalues[1] ) <ftol ) && ( eigenvalues[1] > eigenvalues[2] || fabs ( eigenvalues[1]-eigenvalues[2] ) <ftol );


    return ( check_validity_Q );
}

void mohrcoulomb::ComputePlanePrincipalStressDeriv ( NRmatrix<Doub> &DPSTRS, Doub &epsbarp,NRvector<Int> order ) const
{
}

void mohrcoulomb::ComputePlaneDep ( NRvector<Doub> strain_trial,NRvector<Doub> sigma_proj,NRmatrix<Doub> &eigenvec, Doub &epsbarp,NRvector<Int> order,NRmatrix<Doub> &Dep ) const
{
}

void mohrcoulomb::ComputePlaneTangent ( NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order ) const
{
}


template<class T>
bool mohrcoulomb::ReturnMapLeftEdge ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                                      TComputeSequence &memory, Doub &epsbarnew ) const
{
    sigma_projected = sigma_trial;
    NRvector<T>  eigenvalues = sigma_projected;
    const Doub sinphi = sin ( fPhi );
    const Doub sinpsi = sin ( fPsi );
    const Doub cosphi = cos ( fPhi );
    const Doub sinphi2 = sinphi*sinphi;
    const Doub cosphi2 = 1. - sinphi2;
    NRvector<T> gamma ( 2, 0. ), phi ( 2, 0. ), sigma_bar ( 2, 0. ), ab ( 2, 0. );
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];

    NRvector<T> phival ( 2, 0. );
    NRvector<NRvector<T>> jac ( 2 ), jac_inv ( 2 );
    for ( int i = 0; i < 2; i++ ) {
        jac[i].assign ( 2, 0. );
        jac_inv[i].assign ( 2, 0. );
    }

    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+ ( eigenvalues[0] + eigenvalues[2] ) * T ( sinphi );
    sigma_bar[1] = eigenvalues[1] - eigenvalues[2]+ ( eigenvalues[1] + eigenvalues[2] ) * T ( sinphi );
    T c, H;
    ab[0] = T ( 4. * fmu* ( 1 + sinphi * sinpsi / 3. ) + 4. * fK* sinphi * sinpsi );
    ab[1] = T ( 2. * fmu* ( 1. - sinphi - sinpsi - sinphi * sinpsi / 3. ) + 4. * fK * sinphi * sinpsi );
    T epsbar = T ( fEpsPlasticBar ) + ( gamma[0] + gamma[1] ) * T ( 2. * cosphi );
    PlasticityFunction ( epsbar, c, H );

    phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T ( 2. * cosphi ) * c;
    phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T ( 2. * cosphi ) * c;
    T res = ( fabs ( phival[0] ) + fabs ( phival[1] ) );
    int n_iterations = 300;
    int i;
    bool stop_criterion;
    for ( i = 0; i < n_iterations; i++ ) {

        jac[0][0] = -ab[0] - T ( 4. * cosphi2 ) * H;
        jac[1][0] = -ab[1] - T ( 4. * cosphi2 ) * H;
        jac[0][1] = -ab[1] - T ( 4. * cosphi2 ) * H;
        jac[1][1] = -ab[0] - T ( 4. * cosphi2 ) * H;

        T det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];


        if ( fabs ( det_jac ) <ftol ) {
            std::cerr << "MohrCoulomb:: Singular jacobian." << std::endl;
            DebugStop();
        }

        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;

        gamma[0] -= ( jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1] );
        gamma[1] -= ( jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1] );

        epsbar = T ( fEpsPlasticBar )+ ( gamma[0] + gamma[1] ) * T ( 2. * cosphi );
        PlasticityFunction ( epsbar, c, H );

        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T ( 2. * cosphi ) * c;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T ( 2. * cosphi ) * c;
        res = ( fabs ( phi[0] ) + fabs ( phi[1] ) );

        if ( fabs ( res ) <ftol ) {
            break;
        }
    }

    if ( i == n_iterations ) {
        //DebugStop();
    }

    memory.fGamma[0] = gamma[0];
    memory.fGamma[1] = gamma[1];
    eigenvalues[0] += -T ( 2. * fmu* ( 1 + sinpsi / 3. ) + 2. * fK * sinpsi ) * gamma[0] + T ( ( 4. * fmu / 3. - 2. * fK ) * sinpsi ) * gamma[1];
    eigenvalues[1] += T ( ( 4. * fmu / 3. - fK*2. ) * sinpsi ) * gamma[0] - T ( 2. * fmu* ( 1. + sinpsi / 3. ) + 2. * fK * sinpsi ) * gamma[1];
    eigenvalues[2] += T ( 2. * fmu* ( 1 - sinpsi / 3. ) - 2. * fK * sinpsi ) * ( gamma[0] + gamma[1] );
    sigma_projected = eigenvalues;
    epsbarnew = epsbar;


    bool check_validity_Q = ( eigenvalues[0] > eigenvalues[1] || fabs ( eigenvalues[0]-eigenvalues[1] ) <ftol ) && ( eigenvalues[1] > eigenvalues[2] || fabs ( eigenvalues[1]-eigenvalues[2] ) <ftol );
    return ( check_validity_Q );
    return ( check_validity_Q );
}

void mohrcoulomb::ComputeLeftEdgeTangent ( NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order ) const
{
}

template<class T>
bool mohrcoulomb::ReturnMapRightEdge ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                                       TComputeSequence &memory, Doub &epsbarnew ) const
{
    sigma_projected = sigma_trial;
    NRvector<T> eigenvalues = sigma_projected;
    const Doub sinphi = sin ( fPhi );
    const Doub sinpsi = sin ( fPsi );
    const Doub cosphi = cos ( fPhi );
    const Doub sinphi2 = sinphi*sinphi;
    const Doub cosphi2 = 1. - sinphi2;
    const Doub KV = fK;
    const Doub GV = fmu;

    NRvector<T> gamma ( 2, 0. ), phi ( 2, 0. ), sigma_bar ( 2, 0. ), ab ( 2, 0. );
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    NRvector<T> phival ( 2, 0. );
    NRvector<NRvector<T>> jac ( 2 ), jac_inv ( 2 );
    for ( int i = 0; i < 2; i++ ) {
        jac[i].assign ( 2, 0. );
        jac_inv[i].assign ( 2, 0. );
    }

    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+ ( eigenvalues[0] + eigenvalues[2] ) * T ( sinphi );
    sigma_bar[1] = eigenvalues[0] - eigenvalues[1]+ ( eigenvalues[0] + eigenvalues[1] ) * T ( sinphi );
    T c, H;
    ab[0] = T ( 4. * GV * ( 1 + sinphi * sinpsi / 3. ) + 4. * KV * sinphi * sinpsi );
    ab[1] = T ( 2. * GV * ( 1. + sinphi + sinpsi - sinphi * sinpsi / 3. ) + 4. * KV * sinphi * sinpsi );
    T epsbar = T ( fEpsPlasticBar )+ ( gamma[0] + gamma[1] ) * T ( 2. * cosphi );
    PlasticityFunction ( epsbar, c, H );

    phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T ( 2. * cosphi ) * c;
    phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T ( 2. * cosphi ) * c;

    T res = ( fabs ( phival[0] ) + fabs ( phival[1] ) );
    int n_iterations = 300; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for ( i = 0; i < n_iterations; i++ ) {

        jac[0][0] = -ab[0] - T ( 4. * cosphi2 ) * H;
        jac[1][0] = -ab[1] - T ( 4. * cosphi2 ) * H;
        jac[0][1] = -ab[1] - T ( 4. * cosphi2 ) * H;
        jac[1][1] = -ab[0] - T ( 4. * cosphi2 ) * H;

        T det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

        if ( fabs ( det_jac ) <ftol ) {
            std::cerr << "MohrCoulombPV:: Singular jacobian." << std::endl;
            DebugStop();
        }


        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;

        gamma[0] -= ( jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1] );
        gamma[1] -= ( jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1] );

        epsbar = T ( fEpsPlasticBar )+ ( gamma[0] + gamma[1] ) * T ( 2. * cosphi );
        PlasticityFunction ( epsbar, c, H );

        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - T ( 2. * cosphi ) * c;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - T ( 2. * cosphi ) * c;
        phival[0] = phi[0];
        phival[1] = phi[1];
        res = ( fabs ( phival[0] ) + fabs ( phival[1] ) );


        if ( fabs ( res ) <ftol ) {
            break;
        }
    }

    if ( i == n_iterations ) {
        //DebugStop();
    }


    memory.fGamma[0] = gamma[0];
    memory.fGamma[1] = gamma[1];

    eigenvalues[0] -= T ( 2. * GV * ( 1 + sinpsi / 3. ) + 2. * KV * sinpsi ) * ( gamma[0] + gamma[1] );
    eigenvalues[1] += T ( ( 4. * GV / 3. - KV * 2. ) * sinpsi ) * gamma[0] + T ( 2. * GV * ( 1. - sinpsi / 3. ) - 2. * KV * sinpsi ) * gamma[1];
    eigenvalues[2] += T ( 2. * GV * ( 1 - sinpsi / 3. ) - 2. * KV * sinpsi ) * gamma[0] + T ( ( 4. * GV / 3. - 2. * KV ) * sinpsi ) * gamma[1];
    sigma_projected = eigenvalues;
    epsbarnew = epsbar;

    bool check_validity_Q = ( eigenvalues[0] > eigenvalues[1] || fabs ( eigenvalues[0]-eigenvalues[1] ) <ftol ) && ( eigenvalues[1]  >  eigenvalues[2]  || fabs ( eigenvalues[1]-eigenvalues[2] ) <ftol );
    return ( check_validity_Q );
}

void mohrcoulomb::ComputeRightEdgeTangent ( NRmatrix<Doub> &tang, Doub &epsbarp,  NRvector<Int> order ) const
{
}

template<class T>
bool mohrcoulomb::ReturnMapApex ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                                  TComputeSequence &memory, Doub &epsbarnew ) const
{
    const Doub K = fK;
    const Doub sinpsi = sin ( fPsi );
    const Doub cosphi = cos ( fPhi );
    const Doub cotphi = 1. / tan ( fPhi );
    T ptrnp1 = 0.;
    for ( int i = 0; i < 3; i++ ) {
        ptrnp1 += T ( sigma_trial[i] );
    }
    ptrnp1 /= 3.;
    T DEpsPV = 0.;
    T epsbarnp1 = T ( fEpsPlasticBar );
    T c, H;
    PlasticityFunction ( epsbarnp1, c, H );

    T alpha = cos ( fPhi ) / sin ( fPsi );
    Doub tol = ftol;

    T res = c * cotphi - ptrnp1;
    T pnp1;

    int n_iterations = 300; // @TODO : Define a numeric controls manager object and use it to obtain this information
    int i;
    bool stop_criterion;
    for ( i = 0; i < n_iterations; i++ ) {
        const T jac = H * T ( cosphi * cotphi ) / T ( sinpsi ) + T ( K );
        DEpsPV -= res / jac;

        epsbarnp1 = T ( fEpsPlasticBar ) + T ( alpha ) * DEpsPV;
        pnp1 = ptrnp1 - T ( K ) * DEpsPV;
        PlasticityFunction ( epsbarnp1, c, H );
        res = c * cotphi - pnp1;

        if ( fabs ( res ) <ftol ) {
            break;
        }
    }

    if ( i == n_iterations ) {
        //  DebugStop();
    }

    epsbarnew = epsbarnp1;
    for ( int i = 0; i < 3; i++ ) {
        sigma_projected[i] = pnp1;
    }
    return true;
}


void mohrcoulomb::ComputeApexGradient ( NRmatrix<Doub> & gradient, Doub & eps_bar_p, NRvector<Int>  order ) const
{
}


void mohrcoulomb::ProjectSigma ( const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub>  &gradient, NRvector<Int> order,NRmatrix<Doub> &eigenvec,NRvector<Doub> pstrain )
{



//   cout << "BULK "<<fK <<endl;
    gradient.assign ( 3,3,0. );
    TComputeSequence memory;
    this->SetEpsBar ( k_prev );
    Doub epsbartemp = k_prev; // it will be defined by the correct returnmap

    bool check_validity_Q;

    // Check if we are in the correct sextant
    check_validity_Q = ( sigma_trial[0] > sigma_trial[1] || fabs ( sigma_trial[0]-sigma_trial[1] ) <ftol ) && ( sigma_trial[1] > sigma_trial[2] || fabs ( sigma_trial[1]-sigma_trial[2] ) <ftol );
    if ( !check_validity_Q ) {
        DebugStop();
    }

//22652.243938097836
    Doub phi = PhiPlane<Doub> ( sigma_trial );
    bool elastic_update_Q = fabs ( phi ) <ftol || phi < 0.0;
    if ( elastic_update_Q ) {
        m_type = 0; // Elastic behavior
        memory.fWhichPlane = TComputeSequence::EElastic;
        memory.fGamma.resize ( 0 );
        sigma = sigma_trial;

        NRmatrix<Doub> C = GetElasticMatrix();
        gradient[0][0] = C[0][0];
        gradient[0][1] = C[0][1];
        gradient[0][2] = C[0][5];
        gradient[1][0] = C[1][0];
        gradient[1][1] = C[1][1];
        gradient[1][2] = C[1][5];
        gradient[2][0] = C[5][0];
        gradient[2][1] = C[5][1];
        gradient[2][2] = C[5][5];

        return;
    }

    m_type = 1; // failure behavior
    NRvector<Doub>sigma_projected;
    memory.fGamma.resize ( 1 );
    memory.fGamma[0] = 0.;
    check_validity_Q = this->ReturnMapPlane<Doub> ( sigma_trial, sigma_projected, memory, epsbartemp );
    if ( check_validity_Q ) {
        k_proj = epsbartemp;
        this->SetEpsBar ( k_proj );
        sigma = sigma_projected;
        memory.fWhichPlane = TComputeSequence::EMainPlane;
        ComputePlaneTangent ( gradient, epsbartemp,order );
        ComputePlaneDep ( pstrain,sigma_projected,eigenvec, epsbartemp,order,gradient );
        //cout << "plane"<<endl;

    } else {
        memory.fGamma.resize ( 2 );
        memory.fGamma[0] = 0.;
        memory.fGamma[1] = 0.;
        bool IsEdge = false;

        //cout << "?"<<endl;
        const Doub sinpsi = sin ( fPsi );
        Doub val = ( 1 - sinpsi ) * sigma_trial[0] - 2. * sigma_trial[1] + ( 1 + sinpsi ) * sigma_trial[2];
        if ( val > 0. ) {
            IsEdge = this->ReturnMapRightEdge<Doub> ( sigma_trial, sigma_projected, memory, epsbartemp );
            memory.fWhichPlane = TComputeSequence::ERightEdge;
            //cout << "rigth"<<endl;
            ComputeRightEdgeTangent ( gradient, epsbartemp, order );
            //cout << "rigth"<<endl;

        } else {
            IsEdge = this->ReturnMapLeftEdge<Doub> ( sigma_trial, sigma_projected, memory, epsbartemp );
            memory.fWhichPlane = TComputeSequence::ELeftEdge;
            //cout << "left"<<endl;

            ComputeLeftEdgeTangent ( gradient, epsbartemp,order );
            //cout << "left"<<endl;
        }
        if ( !IsEdge ) {
            m_type = -1; // Tensile behavior
            this->ReturnMapApex ( sigma_trial, sigma_projected, memory, epsbartemp );
            memory.fWhichPlane = TComputeSequence::EApex;
            //cout << "apex"<<endl;


            ComputeApexGradient ( gradient, epsbartemp,order );
        }

        k_proj = epsbartemp;
        this->SetEpsBar ( k_proj );
        sigma = sigma_projected;
    }

    //cout<<"PROJ"<<endl;
    //  sigma_projected.Print();


    if ( fsetconsistentangent==false ) {
        NRmatrix<Doub> C = GetElasticMatrix();
        gradient[0][0] = C[0][0];
        gradient[0][1] = C[0][1];
        gradient[0][2] = C[0][5];
        gradient[1][0] = C[1][0];
        gradient[1][1] = C[1][1];
        gradient[1][2] = C[1][5];
        gradient[2][0] = C[5][0];
        gradient[2][1] = C[5][1];
        gradient[2][2] = C[5][5];
    }

}

void mohrcoulomb::updateatributes ( NRvector<MatDoub> mult )
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
    SetUp ( fyoung,fnu, fc,fPhi, fPsi );
}

bool  mohrcoulomb::ProjectSigma ( NRvector<Doub> &sigma_trial, Doub &k_prev,NRvector<Doub> &sigma,Doub &k_proj, Int & whatphi )
{
    //   cout << "BULK "<<fK <<endl;
    TComputeSequence memory;
    this->SetEpsBar ( k_prev );
    Doub epsbartemp = k_prev; // it will be defined by the correct returnmap

    bool check_validity_Q;

    // Check if we are in the correct sextant
    check_validity_Q = ( sigma_trial[0] > sigma_trial[1] || fabs ( sigma_trial[0]-sigma_trial[1] ) <ftol ) && ( sigma_trial[1] > sigma_trial[2] || fabs ( sigma_trial[1]-sigma_trial[2] ) <ftol );
    if ( !check_validity_Q ) {
        DebugStop();
    }
    //
//22652.243938097836
    Doub phi = PhiPlane<Doub> ( sigma_trial );
    bool elastic_update_Q = fabs ( phi ) <ftol || phi < 0.0;
    if ( elastic_update_Q ) {
        memory.fWhichPlane = TComputeSequence::EElastic;
        memory.fGamma.resize ( 0 );
        sigma = sigma_trial;


        return elastic_update_Q;
    }

    NRvector<Doub>sigma_projected;
    memory.fGamma.resize ( 1 );
    memory.fGamma[0] = 0.;
    check_validity_Q = this->ReturnMapPlane<Doub> ( sigma_trial, sigma_projected, memory, epsbartemp );
    if ( check_validity_Q ) {
        k_proj = epsbartemp;
        this->SetEpsBar ( k_proj );
        sigma = sigma_projected;
        memory.fWhichPlane = TComputeSequence::EMainPlane;
        //cout<<"adsas"<<endl;
        //cout << "plane"<<endl;
        whatphi=0;

    } else {

        memory.fGamma.resize ( 2 );
        memory.fGamma[0] = 0.;
        memory.fGamma[1] = 0.;
        bool IsEdge = false;

        //cout << "?"<<endl;
        const Doub sinpsi = sin ( fPsi );
        Doub val = ( 1 - sinpsi ) * sigma_trial[0] - 2. * sigma_trial[1] + ( 1 + sinpsi ) * sigma_trial[2];
        if ( val > 0. ) {
            IsEdge = this->ReturnMapRightEdge<Doub> ( sigma_trial, sigma_projected, memory, epsbartemp );
            memory.fWhichPlane = TComputeSequence::ERightEdge;
            //cout << "rigth"<<endl;
            whatphi=1;
        } else {
            IsEdge = this->ReturnMapLeftEdge<Doub> ( sigma_trial, sigma_projected, memory, epsbartemp );
            memory.fWhichPlane = TComputeSequence::ELeftEdge;
            //cout << "left"<<endl;
            whatphi=2;
        }
        if ( !IsEdge ) {
            this->ReturnMapApex ( sigma_trial, sigma_projected, memory, epsbartemp );
            memory.fWhichPlane = TComputeSequence::EApex;
            //cout << "apex"<<endl;
            whatphi=3;
        }

        k_proj = epsbartemp;
        this->SetEpsBar ( k_proj );
        sigma = sigma_projected;
    }
    return elastic_update_Q;
}

bool mohrcoulomb::closestpointproj ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain,Int &whatphi )
{
    bool iselastic=true;
    NRtensor<Doub>  epse = epst - epsp;
    NRmatrix<Doub>  C = GetElasticMatrix();

    NRmatrix<Doub>  tempepsemat, stresstrial, Dept;
    epse.FromTensorToNRmatrix ( tempepsemat );

    C.Mult ( tempepsemat, stresstrial );

    NRtensor<Doub>  stresstrialtensor,epsptensor;
    epse.FromNRmatrixToTensor ( stresstrial, stresstrialtensor );





    NRvector<Doub> eigenvaluesinit;
    epst.EigenValue ( eigenvaluesinit );

    NRmatrix<Doub>  pt, vec;
    NRvector<Doub> sigma_trial ( 3,0. ),sigma,gradient ( 3,3 );

    stresstrialtensor.EigenSystem ( pt, vec );


    sigma_trial[0] = pt[0][0];
    sigma_trial[1] =  pt[0][1];
    sigma_trial[2] = pt[0][2];

    Doub k_prev=0.;
    Doub k_proj=0.;
    Int  m_type;
    NRvector<Doub> princstrain ( 3 );


    iselastic=ProjectSigma ( sigma_trial, k_prev, sigma, k_proj,whatphi );

    //sigma.Print();

    NRmatrix<Doub>  fulltensorproj = stressrecosntruction ( sigma, vec );


    NRmatrix<Doub>  projVoigtMat, invCe = GetInverseElasticMatrix(), epsemat;
    fulltensorproj.FromFullToVoigt ( projVoigtMat );

    NRtensor<Doub>  voigtProjTensor;
    voigtProjTensor.FromNRmatrixToTensor ( projVoigtMat, voigtProjTensor );
    projstress = voigtProjTensor;
    //NRmatrix<Doub> Dep;
    //ComputePlaneTangent ( projstress, stresstrialtensor,  Dep );

    //Dep.Print();
    invCe.Mult ( projVoigtMat, epsemat );
    NRmatrix<Doub>  epspmat = tempepsemat;
    epspmat -= epsemat;

    voigtProjTensor.FromNRmatrixToTensor ( epspmat, epsptensor );

    voigtProjTensor.FromNRmatrixToTensor ( epsemat, projstrain );


    return iselastic;
}

void mohrcoulomb::closestpointproj2 ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma )
{
    NRtensor<Doub>  epse = epst - epsp;
    NRmatrix<Doub>  C = GetElasticMatrix();

    NRmatrix<Doub>  tempepsemat, stresstrial, Dept;
    epse.FromTensorToNRmatrix ( tempepsemat );

    C.Mult ( tempepsemat, stresstrial );

    NRtensor<Doub>  stresstrialtensor,epsptensor;
    epse.FromNRmatrixToTensor ( stresstrial, stresstrialtensor );



    NRvector<Doub> eigenvaluesinit;
    epst.EigenValue ( eigenvaluesinit );

    cout <<"sdsdas"<<endl;
    stresstrial.Print();

    NRvector<Int> order ( 3 );

    order[0]=0;
    order[1]=2;
    order[2]=1;
    //cout <<"EIGX"<<endl;
    NRvector<Doub> PSTRS ( 3 );
    Doub TRX=stresstrialtensor.XX()+stresstrialtensor.YY();
    Doub B =sqrt ( pow ( ( stresstrialtensor.XX()-stresstrialtensor.YY() ),2 )+4.*stresstrialtensor.XY() *stresstrialtensor.XY() );
    PSTRS[0]=0.5* ( TRX+B );
    PSTRS[1]=0.5* ( TRX-B );
    PSTRS[2]=2.*fmu* ( epst.ZZ()-epst.I1() /3. )+fK*epst.I1();
    //PSTRS.Print();


    // Identify maximum (PSTRS1) and minimum (PSTRS3) principal stresses
    Doub II=0;
    Doub JJ=0;
    Doub MM;
    Doub PSTRS1=PSTRS[II];
    Doub PSTRS3=PSTRS[JJ];
    for ( Int I = 1; I <= 2; I++ ) {
        if ( PSTRS[I]>=PSTRS1 ) {
            II=I;
            PSTRS1=PSTRS[II];
        }
        if ( PSTRS[I]<PSTRS3 ) {
            JJ=I;
            PSTRS3=PSTRS[JJ];
        }
    }
    if ( II!=0&&JJ!=0 ) {
        MM=0;
    }
    if ( II!=1&&JJ!=1 ) {
        MM=1;
    }
    if ( II!=2&&JJ!=2 ) {
        MM=2;
    }
    Doub PSTRS2=PSTRS[MM];
    order[0]=II;
    order[1]=JJ;
    order[2]=MM;

    /* order[0]=0;
     order[1]=1;
     order[2]=2;*/

    cout<<"val e vec"<<endl;
    //order.Print();


    NRmatrix<Doub>  pt, vec,pt2, vec2;
    NRvector<Doub> sigma_trial ( 3,0. ),sigma,gradient ( 3,3 );

    stresstrialtensor.EigenSystem ( pt, vec );
    epst.EigenSystem ( pt2, vec2 );

    pt2.Print();
    vec2.Print();
    sigma_trial[0] = pt[0][0];
    sigma_trial[1] =  pt[0][1];
    sigma_trial[2] = pt[0][2];

    cout<<"sigma_trial"<<endl;
    sigma_trial.Print();

    Doub k_prev=0.;
    Doub k_proj=0.;
    Int  m_type;
    NRvector<Doub> princstrain ( 3 );
    princstrain[0] = pt2[0][0];
    princstrain[1] =  pt2[0][1];
    princstrain[2] = pt2[0][2];
    ProjectSigma ( sigma_trial, k_prev, sigma, k_proj,m_type, Dep,order,vec,princstrain );
    cout<<"sigma_proj"<<endl;
    sigma.Print();

    NRmatrix<Doub>  fulltensorproj = stressrecosntruction ( sigma, vec );


    NRmatrix<Doub>  projVoigtMat, invCe = GetInverseElasticMatrix(), epsemat;
    fulltensorproj.FromFullToVoigt ( projVoigtMat );

    NRtensor<Doub>  voigtProjTensor;
    voigtProjTensor.FromNRmatrixToTensor ( projVoigtMat, voigtProjTensor );
    projstress = voigtProjTensor;


    //cout << "nvec = " << endl;
    //nvec.Print();

    invCe.Mult ( projVoigtMat, epsemat );
    NRmatrix<Doub>  epspmat = tempepsemat;
    epspmat -= epsemat;

    voigtProjTensor.FromNRmatrixToTensor ( epspmat, epsptensor );

    voigtProjTensor.FromNRmatrixToTensor ( epsemat, projstrain );


}

void mohrcoulomb::closestpointproj ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep2, Doub & projgamma )
{
    NRtensor<Doub>  epse = epst - epsp;
    NRmatrix<Doub>  C = GetElasticMatrix();

    NRmatrix<Doub>  tempepsemat, stresstrial, Dept;
    epse.FromTensorToNRmatrix ( tempepsemat );

    C.Mult ( tempepsemat, stresstrial );

    NRtensor<Doub>  stresstrialtensor,epsptensor;
    epse.FromNRmatrixToTensor ( stresstrial, stresstrialtensor );



    NRvector<Doub> eigenvaluesinit;
    epst.EigenValue ( eigenvaluesinit );


    NRmatrix<Doub>  pt, vec,pt2, vec2;
    NRvector<Doub> sigma_trial ( 3,0. ),sigma,gradient ( 3,3 );

    stresstrialtensor.EigenSystem ( pt, vec );
    epst.EigenSystem ( pt2, vec2 );

    sigma_trial[0] = pt[0][0];
    sigma_trial[1] =  pt[0][1];
    sigma_trial[2] = pt[0][2];

    Doub k_prev=0.;
    Doub k_proj=0.;
    Int  m_type;
    NRvector<Doub> princstrain ( 3 );
    princstrain[0] = pt2[0][0];
    princstrain[1] =  pt2[0][1];
    princstrain[2] = pt2[0][2];
    bool iselastic=ProjectSigma ( sigma_trial, k_prev, sigma, k_proj,m_type );
    //cout << "tipo" << endl;
    //cout << m_type << endl;

    NRmatrix<Doub>  fulltensorproj = stressrecosntruction ( sigma, vec );


    NRmatrix<Doub>  projVoigtMat, invCe = GetInverseElasticMatrix(), epsemat;
    fulltensorproj.FromFullToVoigt ( projVoigtMat );

    NRtensor<Doub>  voigtProjTensor;
    voigtProjTensor.FromNRmatrixToTensor ( projVoigtMat, voigtProjTensor );
    projstress = voigtProjTensor;


    invCe.Mult ( projVoigtMat, epsemat );
    NRmatrix<Doub>  epspmat = tempepsemat;
    epspmat -= epsemat;

    voigtProjTensor.FromNRmatrixToTensor ( epspmat, epsptensor );

    voigtProjTensor.FromNRmatrixToTensor ( epsemat, projstrain );
    Dep2.assign ( 3,3,0. );
    NRmatrix<Doub> Dep ( 6,6,0. );

// 	NRmatrix<Doub> nvec = avec ( stresstrialtensor );
//
// 	NRvector<NRmatrix<Doub>> nvecs = N( stresstrialtensor );
//
// 	//
// 	NRmatrix<Doub> nvec1 = nvecs[0];
// 	NRmatrix<Doub> n1,temp,temp2,tempstresstriela(stresstrial);
// 	tempstresstriela-=projVoigtMat;
// 	C.Mult(nvec,temp);
// 	C.Mult(nvec1,temp2);
//
// 	 Doub gamma3 = tempstresstriela.NRmatrixNorm() / temp.NRmatrixNorm();
// 	 Doub gamma4 = tempstresstriela.NRmatrixNorm() / temp2.NRmatrixNorm();
//
//
// 	cout << "gamma 3333" << endl;
// 	cout << gamma3  << endl;
//
// 	cout << "gamma 4444" << endl;
// 	cout << gamma4  << endl;
// 	cout << gamma2  << endl;
	
	
	    

    if ( iselastic==true ) {
        //cout << "elastic"  << endl;
        Dep2[0][0] = C[0][0];
        Dep2[0][1] = C[0][1];
        Dep2[0][2] = C[0][5];
        Dep2[1][0] = C[1][0];
        Dep2[1][1] = C[1][1];
        Dep2[1][2] = C[1][5];
        Dep2[2][0] = C[5][0];
        Dep2[2][1] = C[5][1];
        Dep2[2][2] = C[5][5];

        return;
    }

    if ( m_type==1 || m_type==2 ) {
        //cout << "corner"  << endl;
        ComputeConsistentEdgeTangent ( projstress, stresstrialtensor, Dep,m_type );
        Dep2[0][0] = Dep[0][0];
        Dep2[0][1] = Dep[0][1];
        Dep2[0][2] = Dep[0][5];
        Dep2[1][0] = Dep[1][0];
        Dep2[1][1] = Dep[1][1];
        Dep2[1][2] = Dep[1][5];
        Dep2[2][0] = Dep[5][0];
        Dep2[2][1] = Dep[5][1];
        Dep2[2][2] = Dep[5][5];
	//closestpointproj3 (  epst,  epsp,projstress,projstrain,Dep2, projgamma );
        return;
    } else if ( m_type==0 ) {
        //cout << "plane"  << endl;
        ComputeConsistentPlaneTangent ( projstress, stresstrialtensor, Dep );
        Dep2[0][0] = Dep[0][0];
        Dep2[0][1] = Dep[0][1];
        Dep2[0][2] = Dep[0][5];
        Dep2[1][0] = Dep[1][0];
        Dep2[1][1] = Dep[1][1];
        Dep2[1][2] = Dep[1][5];
        Dep2[2][0] = Dep[5][0];
        Dep2[2][1] = Dep[5][1];
        Dep2[2][2] = Dep[5][5];
		//closestpointproj3 (  epst,  epsp,projstress,projstrain,Dep2, projgamma );
        return;
    } else if ( m_type==3 ) { //Apex
        //cout << "apex"  << endl;
		//closestpointproj3 (  epst,  epsp,projstress,projstrain,Dep2, projgamma );
        Dep2[0][0] =1.e-12;
        Dep2[0][1] = 1.e-12;
        Dep2[0][2] = 1.e-12;
        Dep2[1][0] = 1.e-12;
        Dep2[1][1] = 1.e-12;
        Dep2[1][2] = 1.e-12;
        Dep2[2][0] =1.e-12;
        Dep2[2][1] = 1.e-12;
        Dep2[2][2] = 1.e-12;
//closestpointproj3 (  epst,  epsp,projstress,projstrain,Dep2, projgamma );
// 				Dep2[0][0] = C[0][0];Dep2[0][1] = C[0][1];Dep2[0][2] = C[0][5];
// 		Dep2[1][0] = C[1][0];Dep2[1][1] = C[1][1];Dep2[1][2] = C[1][5];
// 		Dep2[2][0] = C[5][0];Dep2[2][1] = C[5][1];Dep2[2][2] = C[5][5];
        return;
    }



}

void mohrcoulomb::ComputeConsistentPlaneTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep )
{

    Doub J2=trialstress.J2();
    Doub J3=trialstress.J3();
    Doub val = -3*sqrt ( 3. ) *J3/ ( 2.*pow ( J2,1.5 ) );
    if ( val>1. ) {
        val=1.;
    }
    if ( val<-1. ) {
        val=-1.;
    }
    ftheta= 1/3.*asin ( val );

    fa = cos ( ftheta )-1./sqrt ( 3. ) *sin ( ftheta ) *sin ( fPhi );

    fda=- ( ( cos ( ftheta ) *sin ( fPhi ) ) /sqrt ( 3 ) ) - sin ( ftheta );

    fd2a = -cos ( ftheta ) + ( sin ( fPhi ) *sin ( ftheta ) ) /sqrt ( 3 );


    //Compute Deltalambda
    NRmatrix<Doub> E=GetElasticMatrix();
    NRmatrix<Doub> atemp,denom;
    Doub phiinvars=PhiInvars ( trialstress );
	//Doub phiinvars=Yields(trialstress)[0];
    NRmatrix<Doub> a,at;
	a=avec ( trialstress );
	//a =N(trialstress)[0];
    a.Transpose ( at );
    at.Mult ( E,atemp );
    atemp.Mult ( a,denom );
    Doub Deltalambda = phiinvars/denom[0][0];

    //cout << "Deltalambda"<< Deltalambda << endl;
    NRmatrix<Doub> dadsig = dAdsig ( trialstress ),Etmodular ( E ),InverseTemp,Q;
    Q.IdentityMatrix ( 6 );

    NRmatrix<Doub> Ea,atEt,atEa,Et,atE,num;
    E.Transpose ( Et );
    E.Mult ( a,Ea );
    at.Mult ( Et,atEt );
    at.Mult ( E,atE );
    atE.Mult ( a,atEa );
    Ea.Mult ( atEt,num );
    num*=1/atEa[0][0];
    Etmodular-=num;

    NRmatrix<Doub> temp;

    dadsig.Mult ( E,temp );
    temp*=Deltalambda;

    Q-=temp;

    Etmodular.Mult ( Q,Dep );
}

void mohrcoulomb::ComputeConsistentEdgeTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep, Int &m_type )
{
    Doub J2=trialstress.J2();
    Doub J3=trialstress.J3();
    Doub val = -3*sqrt ( 3. ) *J3/ ( 2.*pow ( J2,1.5 ) );
    if ( val>1. ) {
        val=1.;
    }
    if ( val<-1. ) {
        val=-1.;
    }
    ftheta= 1/3.*asin ( val );

    fa = cos ( ftheta )-1./sqrt ( 3. ) *sin ( ftheta ) *sin ( fPhi );

    fda=- ( ( cos ( ftheta ) *sin ( fPhi ) ) /sqrt ( 3 ) ) - sin ( ftheta );

    fd2a = -cos ( ftheta ) + ( sin ( fPhi ) *sin ( ftheta ) ) /sqrt ( 3 );


    //cout << "theta "<<endl;
    //cout <<ftheta*180/M_PI<<endl;
    // NRvector<NRmatrix<Doub>> nvecs = N ( trialstress );
    // NRvector<Doub> phis = Yields ( trialstress );
    NRmatrix<Doub> E = GetElasticMatrix();
    NRmatrix<Doub> a,b,at,bt;

    a=avec ( trialstress );

    Doub ca,cb,cd,faa,f2a;
    faa=PhiInvars ( trialstress );
    //dadsig para thetax
    //a=nvecs[0];
    //faa=phis[0];

    NRmatrix<Doub> dbdsig,dadsig;
    NRmatrix<Doub> Ea,Eb,atE,btE,tempa,tempb,tempd;
    dadsig = dAdsig ( trialstress );
    if ( m_type==1 ) { //rigth f1 e f2
        // b=nvecs[1];
        //   f2a=phis[1];
        //dbdsig para theta30
        //ftheta = 30*M_PI/180.;
        fa= ( cos ( ftheta ) * ( 1 + sin ( fPhi ) ) ) /2. + ( ( -3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );
        fda= ( cos ( ftheta ) * ( -3 + sin ( fPhi ) ) ) / ( 2.*sqrt ( 3 ) ) - ( ( 1 + sin ( fPhi ) ) *sin ( ftheta ) ) /2.;
        fd2a=- ( cos ( ftheta ) * ( 1 + sin ( fPhi ) ) ) /2. - ( ( -3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );

    } else { //left f1 e f3
        //b=nvecs[2];
        //f2a=phis[2];
        //ftheta = -29.99*M_PI/180.;
        fa= ( cos ( ftheta ) * ( 1 - sin ( fPhi ) ) ) /2. + ( ( 3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );
        fda= ( cos ( ftheta ) * ( 3 + sin ( fPhi ) ) ) / ( 2.*sqrt ( 3 ) ) - ( ( 1 - sin ( fPhi ) ) *sin ( ftheta ) ) /2.;
        fd2a=- ( cos ( ftheta ) * ( 1 - sin ( fPhi ) ) ) /2. - ( ( 3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );
    }

    dbdsig = dAdsig ( trialstress );
    NRmatrix<Doub> symcheck;

    b=avec ( trialstress );
    f2a=PhiInvars ( trialstress );

    a.Transpose ( at );
    b.Transpose ( bt );

    at.Mult ( E,atE );
    bt.Mult ( E,btE );

    atE.Mult ( a,tempa );
    btE.Mult ( b,tempb );
    atE.Mult ( b,tempd );
    Doub q = ( tempa[0][0]*tempb[0][0]-pow ( tempd[0][0],2. ) );
    Doub DeltalambdaA = ( tempb[0][0]*faa-tempd[0][0]*f2a ) / q;
    Doub DeltalambdaB = ( tempa[0][0]*f2a-tempd[0][0]*faa ) / q;

    E.Mult ( a,Ea );
    E.Mult ( b,Eb );

    NRmatrix<Doub>temp1,temp2,temp3,temp4,Et2 ( E );

    Ea.Mult ( atE,temp1 );
    temp1*=tempb[0][0];

    Ea.Mult ( btE,temp2 );
    temp2*=tempd[0][0];

    Eb.Mult ( atE,temp3 );
    temp3*=tempd[0][0];

    Eb.Mult ( btE,temp4 );
    temp4*=tempa[0][0];

    temp1-=temp2;
    temp1-=temp3;
    temp1+=temp4;
    temp1*=1./q;
    Et2-=temp1;
    //Dep=Et2;
    //dadsig.Print();
    //dbdsig.Print();
    //DebugStop();

    NRmatrix<Doub> T,dadsigE,dbdsigE;
    T.IdentityMatrix ( 6 );


    dadsig.Mult ( E,dadsigE );
    //dadsigE.Mult(a,symcheck);
    //cout << "symcheck" << symcheck[0][0] <<endl;

    dbdsig.Mult ( E,dbdsigE );
    //dbdsigE.Mult(b,symcheck);
    //cout << "symcheck B" << symcheck[0][0] <<endl;


    dadsigE*=DeltalambdaA;
    dbdsigE*=DeltalambdaB;

    T-=dadsigE;
    T-=dbdsigE;
    //dbdsig.Print();
    Et2.Mult ( T,Dep );

}

void mohrcoulomb::ComputeConsistentNumericalTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep, Int &m_type )
{

}
