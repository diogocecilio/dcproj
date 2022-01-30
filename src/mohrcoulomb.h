#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <random>
#include "ludcmp.h"
#include <Eigen/Dense>
class mohrcoulomb
{
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
    Doub ftol=1.e-12;
    bool fsetconsistentangent=true;
    Doub ftheta;
    Doub fa,fda,fd2a;


protected:
    Doub fEpsPlasticBar;

public:

    /// structure which contains the decision tree of the return map
    // we can only expect a consistent tangent matrix if the decision tree remains the same

    struct TComputeSequence {

        TComputeSequence() : fWhichPlane ( ENoPlane ), fGamma ( 0 )
        {

        }

        TComputeSequence ( const TComputeSequence &copy ) : fWhichPlane ( copy.fWhichPlane ), fGamma ( copy.fGamma )
        {

        }

        TComputeSequence &operator= ( const TComputeSequence &copy )
        {
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
    mohrcoulomb ( Doub young, Doub nu, Doub c, Doub Phi,Doub Psi );
    mohrcoulomb ( const mohrcoulomb &cp );

    NRmatrix<Doub>  stressrecosntruction ( NRvector<Doub>  val, NRmatrix<Doub>  vec )
    {
        NRmatrix<Doub>  sol ( 3, 3, 0. );
        for ( Int i = 0; i < 3; i++ ) {
            MatDoub colvec ( 3, 1, 0. ), colvect, temp;
            for ( Int j = 0; j < 3; j++ ) {
                colvec[j][0] = vec[j][i];
            }

            colvec.Transpose ( colvect );
            colvec.Mult ( colvect, temp );
            temp *= val[i];
            sol += temp;
        }
        //sol.Print();
        return sol;
    }

    NRmatrix<Doub>  HW ( NRmatrix<Doub>  sig )
    {
        MatDoub sol ( 3, 1 );
        Doub xi, rho, beta;
        xi = sig[0][0];
        rho = sig[1][0];
        beta = sig[2][0];
        sol[0][0] = xi / sqrt ( 3. ) + sqrt ( 2. / 3. ) *rho*cos ( beta );
        sol[1][0] = xi / sqrt ( 3. ) + sqrt ( 2. / 3. ) *rho*cos ( beta - 2. *  M_PI / 3. );
        sol[2][0] = xi / sqrt ( 3. ) + sqrt ( 2. / 3. ) *rho*cos ( beta + 2. *  M_PI / 3. );
        return sol;
    }


    //void SetUp ( Doub Phi, Doub Psi, Doub c,Doub young, Doub nu )
    void SetUp ( Doub young, Doub nu,Doub c,Doub Phi, Doub Psi )
    {
        fPhi = Phi;
        fPsi = Psi;
        fc = c;
        fyoung = young;
        fnu = nu;
        fmu = fyoung/ ( 2. * ( 1. + fnu ) );
        flambda = fyoung*fnu/ ( ( 1. + fnu ) * ( 1. - 2.* fnu ) );
        //cout << "YOUNG = " << fyoung << endl;
        //cout << "fnu" << fnu << endl;
        fK = fyoung / ( 3.* ( 1. - 2.* fnu ) );
    }

    mohrcoulomb & operator= ( const mohrcoulomb &cp );

    void SetEpsBar ( Doub &epsbar )
    {
        fEpsPlasticBar = epsbar;
    }


    void SetElasticResponse ( Doub young, Doub nu )
    {
        fyoung = young;
        fnu = nu;
    }

    bool  ProjectSigma ( NRvector<Doub> &sigma_trial, Doub &k_prev,NRvector<Doub> &sigma,Doub &k_proj,Int & whatphi );

    Doub  InitialDamage ( const NRvector<Doub> &stress_p ) const;

    template <class T>
    void PlasticityFunction ( const T epsp, T &c, T &H ) const;

    template<class T>
    NRvector<T> SigmaElastPV ( const NRvector<T> &deform ) const;

    template<class T>
    T PhiPlane ( const NRvector<T> &sigma ) const;


    template<class T>
    bool ReturnMapPlane ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                          TComputeSequence &memory, Doub &epsbarnew ) const;

    void ComputePlaneTangent ( NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order ) const;

    template<class T>
    bool ReturnMapLeftEdge ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                             TComputeSequence &memory, Doub &epsbarnew ) const;

    void ComputeLeftEdgeTangent ( NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order ) const;

    template<class T>
    bool ReturnMapRightEdge ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                              TComputeSequence &memory, Doub &epsbarnew ) const;

    void ComputeRightEdgeTangent ( NRmatrix<Doub> &tang, Doub &epsbarp, NRvector<Int> order ) const;

    template<class T>
    bool ReturnMapApex ( const NRvector<T> &sigma_trial, NRvector<T> &sigma_projected,
                         TComputeSequence &memory, Doub &epsbarnew ) const;


    void ComputeApexGradient ( NRmatrix<Doub> & gradient, Doub & eps_bar_p, NRvector<Int> order ) const;


    void ProjectSigma ( const NRvector<Doub> & sigma_trial, Doub k_prev, NRvector<Doub> & sigma, Doub &k_proj, Int & m_type, NRmatrix<Doub>  &gradient,NRvector<Int> order,NRmatrix<Doub> &eigenvec,NRvector<Doub> pstrain );


    void Phi ( NRvector<Doub> sig_vec, Doub alpha, NRvector<Doub> &phi ) const;

    Doub Phi()
    {
        return fPhi;
    }


    void ComputePlanePrincipalStressDeriv ( NRmatrix<Doub> &DPSTRS, Doub &epsbarp,NRvector<Int> order ) const;

    void ComputePlaneDep ( NRvector<Doub> strain_trial,NRvector<Doub> sigma_proj,NRmatrix<Doub> &eigenvec, Doub &epsbarp,NRvector<Int> order,NRmatrix<Doub> &Dep ) const;

    Doub PhiInvars ( NRtensor<Doub> tensor )
    {
        Doub I1=tensor.I1();
        Doub J2=tensor.J2();
        Doub thetaval = theta ( tensor );
        Doub atheta = fa;//A(tensor,thetaval);
        Doub thetaint =theta ( tensor );
        Doub s1=2.*sqrt ( J2 ) /sqrt ( 3. ) * ( sin ( thetaint+2.*M_PI/3. ) )+I1/3.;
        //cout << "s1= "<< s1<< endl;
        return 1./3. *I1*sin ( fPhi )+sqrt ( J2 ) *atheta-fc*cos ( fPhi );
    }
    void ComputeConsistentPlaneTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep );
    void ComputeConsistentEdgeTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep, Int &m_type );
    void ComputeConsistentNumericalTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep, Int &m_type );

    void ComputePlaneTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep )
    {

        NRmatrix<Doub>  pt, vec;
        NRvector<Doub> sigma_trial ( 3,0. ),sigma,gradient ( 3,3 ),phiv;
        cout << "trialstress = " << endl;
        trialstress.Print();

        trialstress.EigenSystem ( pt, vec );
        sigma_trial[0] = pt[0][0];
        sigma_trial[1] =  pt[0][1];
        sigma_trial[2] = pt[0][2];
        Doub alpha=0.;
        Doub phi ;

        phi= 1./2.* ( sigma_trial[0]-sigma_trial[2] )+1./2.* ( sigma_trial[0]+sigma_trial[2] ) *sin ( fPhi )-fc*cos ( fPhi );
        //phi= (sigma_trial[0]-sigma_trial[2])+(sigma_trial[0]+sigma_trial[2])*sin(fPhi)-2*fc*cos(fPhi);
        cout << "phi principal = "<< phi << endl;
        cout << "pttt = " << endl;
        pt.Print();
        Doub phiinvars=PhiInvars ( trialstress );
        cout << "phiinvars = "<< phiinvars << endl;
        NRmatrix<Doub> a=avec ( projstress ),at;
        a.Transpose ( at );
        NRmatrix<Doub> E=GetElasticMatrix(),tempe,tempproj,projstress2;
        NRmatrix<Doub> atemp,denom;
        at.Mult ( E,atemp );
        atemp.Mult ( a,denom );

        Doub deltalambda = phiinvars/denom[0][0];

        cout << "gamma Crisfield" << endl;
        cout << deltalambda  << endl;

        NRmatrix<Doub> dadsig = dAdsig ( projstress ),R,temp,temp0,InverseTemp;

        E.Mult ( dadsig,temp0 );

        temp0*=deltalambda;

        temp.IdentityMatrix ( 6 );

        temp+=temp0;

        temp.ComputeInverse ( InverseTemp );

        InverseTemp.Mult ( E,R );

        NRmatrix<Doub> Ra,atRt,atRa,Rt,atR,num;
        R.Transpose ( Rt );
        R.Mult ( a,Ra );
        at.Mult ( Rt,atRt );
        at.Mult ( R,atR );
        atR.Mult ( a,atRa );
        Ra.Mult ( atRt,num );
        num*=1/atRa[0][0];
        R-=num;
        Dep=R;
    }


    void ComputeCornerTangent ( NRtensor<Doub> & projstress, NRtensor<Doub> & trialstress,NRmatrix<Doub> & Dep )
    {

        NRmatrix<Doub>  pt, vec;
        NRvector<Doub> sigma_trial ( 3,0. ),sigma,gradient ( 3,3 ),phiv;

        trialstress.EigenSystem ( pt, vec );
        sigma_trial[0] = pt[0][0];
        sigma_trial[1] =  pt[0][1];
        sigma_trial[2] = pt[0][2];

        Doub phiinvars=PhiInvars ( trialstress );
        NRmatrix<Doub> a=avec ( projstress ),at;
        a.Transpose ( at );
        NRmatrix<Doub> E=GetElasticMatrix(),tempe,tempproj,projstress2;
        NRmatrix<Doub> atemp,denom;
        at.Mult ( E,atemp );
        atemp.Mult ( a,denom );

        Doub deltalambda = phiinvars/denom[0][0];


        NRmatrix<Doub> dadsig = dAdsig ( projstress ),R,temp,temp0,InverseTemp;


        E.Mult ( dadsig,temp0 );

        temp0*=deltalambda;

        temp.IdentityMatrix ( 6 );

        temp+=temp0;

        temp.ComputeInverse ( InverseTemp );

        InverseTemp.Mult ( E,R );

        NRmatrix<Doub> Ra,atRt,atRa,Rt,atR,num;
        R.Transpose ( Rt );
        R.Mult ( a,Ra );
        at.Mult ( Rt,atRt );
        at.Mult ( R,atR );
        atR.Mult ( a,atRa );
        Ra.Mult ( atRt,num );
        num*=1/atRa[0][0];
        R-=num;
        Dep=R;


    }


    Doub theta ( NRtensor<Doub> tensor )
    {
        Doub J2=tensor.J2();
        Doub J3=tensor.J3();
        Doub val = -3*sqrt ( 3. ) *J3/ ( 2.*pow ( J2,1.5 ) );
        if ( val>1. ) {
            val=1.;
        }
        if ( val<-1. ) {
            val=-1.;
        }
        return 1/3.*asin ( val );
    }
    Doub C1()
    {
        //return sin ( fPhi ) /3.;
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return sin ( fPhi ) /3.;
        } else {
            return 0;
        }

    }
    Doub C2 ( NRtensor<Doub> tensor )
    {
        Doub thetaval = ftheta;
        Doub J2=tensor.J2();
        Doub da=fda;//dA(tensor,thetaval);

        Doub a=fa;//A(tensor,thetaval);
        //return 1/2. * pow(J2,-0.5)*(a-tan(3*thetaval)*da);
        //cout << "thetaval = "<< thetaval << endl;
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return 1/2. * pow ( J2,-0.5 ) * ( a-tan ( 3*thetaval ) *da );
        } else {
            return 1/2. * pow ( J2,-0.5 ) *a;
        }
    }
    Doub C3 ( NRtensor<Doub> tensor )
    {
        Doub thetaval = ftheta;
        Doub J2=tensor.J2();
        Doub da=fda;//dA(tensor,thetaval);
        Doub a=fa;//A(tensor,thetaval);
        //return -sqrt(3.)*da/(2.*J2*cos(3*thetaval));
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return -sqrt ( 3. ) *da/ ( 2.*J2*cos ( 3*thetaval ) );
        } else {
            return 0;
        }
    }
    Doub C23 ( NRtensor<Doub> tensor )
    {
        Doub thetaval = ftheta;
        Doub J2=tensor.J2();
        Doub da=fda;//dA(tensor,thetaval);
        Doub a=fa;//A(tensor,thetaval);

        Doub c4 = C4 ( tensor );
        return ( 1./2.*tan ( 3.*thetaval ) *c4+da ) *sqrt ( 3. ) / ( 2*J2*J2*cos ( 3.*thetaval ) );
    }
    Doub C22 ( NRtensor<Doub> tensor )
    {
        Doub thetaval = ftheta;
        Doub J2=tensor.J2();
        Doub da=fda;//dA(tensor,thetaval);

        Doub a=fa;//A(tensor,thetaval);
        Doub c4 = C4 ( tensor );

        return - ( a-pow ( tan ( 3*thetaval ),2 ) *c4-3.*tan ( 3.*thetaval ) *da ) / ( 4.*pow ( J2,3./2. ) );
    }
    Doub C33 ( NRtensor<Doub> tensor )
    {
        Doub J2=tensor.J2();
        Doub thetaval = ftheta;
        Doub c4 = C4 ( tensor );
        return 3.*c4/ ( 4.*pow ( J2,5/2. ) *cos ( 3.*thetaval ) *cos ( 3.*thetaval ) );
    }
    Doub C4 ( NRtensor<Doub> tensor )
    {
        Doub thetaval = ftheta;
        Doub da=fda;//dA(tensor,thetaval);
        Doub d2a=fd2a;//d2A(tensor,thetaval);
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return d2a+3.*tan ( 3.*thetaval ) *da;
        } else {
            return d2a;
        }
    }



    NRmatrix<Doub>  d2J2d2sig()
    {
        NRmatrix<Doub>  P ( 6, 6 );
        P[0][0] = 2. / 3.;
        P[0][1] = -1. / 3.;
        P[0][2] = -1. / 3.;
        P[0][3] = 0.;
        P[0][4] = 0.;
        P[0][5] = 0.;
        P[1][0] = -1. / 3.;
        P[1][1] = 2. / 3.;
        P[1][2] = -1. / 3.;
        P[1][3] = 0.;
        P[1][4] = 0.;
        P[1][5] = 0.;
        P[2][0] = -1. / 3.;
        P[2][1] = -1. / 3.;
        P[2][2] = 2. / 3.;
        P[2][3] = 0.;
        P[2][4] = 0.;
        P[2][5] = 0.;
        P[3][0] = 0.;
        P[3][1] = 0.;
        P[3][2] = 0.;
        P[3][3] = 2.;
        P[3][4] = 0.;
        P[3][5] = 0.;
        P[4][0] = 0.;
        P[4][1] = 0.;
        P[4][2] = 0.;
        P[4][3] = 0.;
        P[4][4] = 2.;
        P[4][5] = 0.;
        P[5][0] = 0.;
        P[5][1] = 0.;
        P[5][2] = 0.;
        P[5][3] = 0.;
        P[5][4] = 0.;
        P[5][5] = 2.;

        return P;
    }

    NRmatrix<Doub>  d2J3d2sig ( NRtensor<Doub> tensor )
    {

        NRmatrix<Doub>  P ( 6, 6 );
        Doub sxx =tensor.XX();
        Doub syy=tensor.YY();
        Doub szz=tensor.ZZ();
        Doub sxy=tensor.XY();
        Doub syz=tensor.YZ();
        Doub sxz=tensor.XZ();

        P[0][0]= ( 4*sxx - 2* ( syy + szz ) ) /9.;
        P[0][1]= ( -2*sxx - 2*syy + 4*szz ) /9.;
        P[0][2]= ( -2*sxx + 4*syy - 2*szz ) /9.;
        P[0][3]= ( 2*sxz ) /3.,P[0][4]= ( -4*syz ) /3.;
        P[0][5]= ( 2*sxy ) /3.;

        P[1][0]= ( -2*sxx - 2*syy + 4*szz ) /9.;
        P[1][1]= ( -2*sxx + 4*syy - 2*szz ) /9.;
        P[1][2]= ( 4*sxx - 2*syy - 2*szz ) /9.;
        P[1][3]= ( -4*sxz ) /3.,P[1][4]= ( 2*syz ) /3.;
        P[1][5]= ( 2*sxy ) /3.;

        P[2][0]= ( -2*sxx + 4*syy - 2*szz ) /9.;
        P[2][1]= ( 4*sxx - 2*syy - 2*szz ) /9.;
        P[2][2]= ( -2*sxx - 2*syy + 4*szz ) /9.;
        P[2][3]= ( 2*sxz ) /3.,P[2][4]= ( 2*syz ) /3.;
        P[2][5]= ( -4*sxy ) /3.;

        P[3][0]= ( 2*sxz ) /3.;
        P[3][1]= ( -4*sxz ) /3.;
        P[3][2]= ( 2*sxz ) /3.;
        P[3][3]= ( 2* ( sxx - 2*syy + szz ) ) /3.,P[3][4]=2*sxy;
        P[3][5]=2*syz;

        P[4][0]= ( -4*syz ) /3.;
        P[4][1]= ( 2*syz ) /3.;
        P[4][2]= ( 2*syz ) /3.;
        P[4][3]=2*sxy,P[4][4]= ( 2* ( -2*sxx + syy + szz ) ) /3.;
        P[4][5]=2*sxz;

        P[5][0]= ( 2*sxy ) /3.;
        P[5][1]= ( 2*sxy ) /3.;
        P[5][2]= ( -4*sxy ) /3.;
        P[5][3]=2*syz,P[5][4]=2*sxz;
        P[5][5]= ( 2* ( sxx + syy - 2*szz ) ) /3.;

        return P;
    }

    NRmatrix<Doub> dAdsig ( NRtensor<Doub> tensor )
    {
        //Eq. 14 (Crisfield Eng. Comput., 1987)
        Doub c1=C1(),c2=C2 ( tensor ),c3=C3 ( tensor ),c4=C4 ( tensor ),c23=C23 ( tensor ),c22=C22 ( tensor ),c33=C33 ( tensor );
        Doub c32 =c23;
        NRmatrix<Doub> da2dsig = d2J2d2sig();
        NRmatrix<Doub> da3dsig = d2J3d2sig ( tensor );
        NRmatrix<Doub> dadsig ( 6,6 ),a2 ( 6,1 ),a3 ( 6,1 ),a2t,a3t,temp1,temp2,temp3,temp4;
        //cout << "dj2"<<endl;
        //tensor.dJ2().Print();
        tensor.dJ2().FromTensorToNRmatrix ( a2 );
        //cout << "a2"<<endl;
        //a2.Print();
        tensor.dJ3().FromTensorToNRmatrix ( a3 );
        a2.Transpose ( a2t );
        a3.Transpose ( a3t ); //tranposto é deitado
        da2dsig*=c2;
        da3dsig*=c3;

        a2.Mult ( a2t,temp1 );
        temp1*=c22;

        a2.Mult ( a3t,temp2 );
        temp2*=c23;

        a3.Mult ( a2t,temp3 );
        temp3*=c32;

        a3.Mult ( a3t,temp4 );
        temp4*=c33;

        dadsig=da2dsig;
        dadsig+=da3dsig;
        dadsig+=temp1;
        dadsig+=temp2;
        dadsig+=temp3;
        dadsig+=temp4;

        return dadsig;
    }

    NRmatrix<Doub> avec ( NRtensor<Doub> tensor )
    {
        //Eq. 14 (Crisfield Eng. Comput., 1987)
        Doub c1=C1(),c2=C2 ( tensor ),c3=C3 ( tensor );
        NRmatrix<Doub> a1 ( 6,1 ),a2 ( 6,1 ),a3 ( 6,1 ),a;
        tensor.dI1().FromTensorToNRmatrix ( a1 );
        tensor.dJ2().FromTensorToNRmatrix ( a2 );
        tensor.dJ3().FromTensorToNRmatrix ( a3 );

        a1*=c1;
        a2*=c2;
        a3*=c3;
        a=a1;
        a+=a2;
        a+=a3;

        return a;
    }

// 	NRtensor<Doub> dA(NRtensor<Doub> tensor)
// 	{
// 		Doub j2=tensor.J2();
// 		Doub j3=tensor.J3();
// 		NRtensor<Doub> dj2=tensor.dJ2();
// 		NRtensor<Doub> dj3=tensor.dJ3();
// 		Doub thetaval = theta(tensor);
// 		Doub tempval=(cos(thetaval)*sin(fPhi) + sqrt(3)*sin(thetaval))/(2.*pow(j2,2.5)*sqrt(4 - (27*pow(j3,2))/pow(j2,3)));
// 		dj2*=(-3.*j3);
// 		dj3*=(2.*j2);
// 		NRtensor<Doub> result(dj2);
// 		result+=dj3;
// 		result*=tempval;
// 		return result;
//
// 	}


    NRvector<NRmatrix<Doub>> N ( NRtensor<Doub> tensor )
    {
        NRvector<NRmatrix<Doub>> out ( 3 );
        NRmatrix<Doub> n1 ( 6,1,0. ),n2 ( 6,1,0. ),n3 ( 6,1,0. );
        NRmatrix<Doub> pt,vec;
        tensor.EigenSystem ( pt,vec );

        Doub e1x=vec[0][0],e1y=vec[1][0],e1z=vec[2][0];
        Doub e3x=vec[0][2],e3y=vec[1][2],e3z=vec[2][2];

        n1[0][0]= ( pow ( e1x,2 ) - pow ( e3x,2 ) ) /2. + ( ( pow ( e1x,2 ) + pow ( e3x,2 ) ) *sin ( fPhi ) ) /2.;
        n1[1][0]= ( pow ( e1y,2 ) - pow ( e3y,2 ) ) /2. + ( ( pow ( e1y,2 ) + pow ( e3y,2 ) ) *sin ( fPhi ) ) /2.;
        n1[2][0]= ( pow ( e1z,2 ) - pow ( e3z,2 ) ) /2. + ( ( pow ( e1z,2 ) + pow ( e3z,2 ) ) *sin ( fPhi ) ) /2.;
        n1[3][0]= ( e1x*e1z - e3x*e3z ) /2. + ( ( e1x*e1z + e3x*e3z ) *sin ( fPhi ) ) /2.;
        n1[4][0]= ( e1y*e1z - e3y*e3z ) /2. + ( ( e1y*e1z + e3y*e3z ) *sin ( fPhi ) ) /2.;
        n1[5][0]= ( e1x*e1y - e3x*e3y ) /2. + ( ( e1x*e1y + e3x*e3y ) *sin ( fPhi ) ) /2.;

        e1x=vec[0][0],e1y=vec[1][0],e1z=vec[2][0];
        e3x=vec[0][1],e3y=vec[1][1],e3z=vec[2][1];

        n2[0][0]= ( pow ( e1x,2 ) - pow ( e3x,2 ) ) /2. + ( ( pow ( e1x,2 ) + pow ( e3x,2 ) ) *sin ( fPhi ) ) /2.;
        n2[1][0]= ( pow ( e1y,2 ) - pow ( e3y,2 ) ) /2. + ( ( pow ( e1y,2 ) + pow ( e3y,2 ) ) *sin ( fPhi ) ) /2.;
        n2[2][0]= ( pow ( e1z,2 ) - pow ( e3z,2 ) ) /2. + ( ( pow ( e1z,2 ) + pow ( e3z,2 ) ) *sin ( fPhi ) ) /2.;
        n2[3][0]= ( e1x*e1z - e3x*e3z ) /2. + ( ( e1x*e1z + e3x*e3z ) *sin ( fPhi ) ) /2.;
        n2[4][0]= ( e1y*e1z - e3y*e3z ) /2. + ( ( e1y*e1z + e3y*e3z ) *sin ( fPhi ) ) /2.;
        n2[5][0]= ( e1x*e1y - e3x*e3y ) /2. + ( ( e1x*e1y + e3x*e3y ) *sin ( fPhi ) ) /2.;

        e1x=vec[0][1],e1y=vec[1][1],e1z=vec[2][1];
        e3x=vec[0][2],e3y=vec[1][2],e3z=vec[2][2];

        n3[0][0]= ( pow ( e1x,2 ) - pow ( e3x,2 ) ) /2. + ( ( pow ( e1x,2 ) + pow ( e3x,2 ) ) *sin ( fPhi ) ) /2.;
        n3[1][0]= ( pow ( e1y,2 ) - pow ( e3y,2 ) ) /2. + ( ( pow ( e1y,2 ) + pow ( e3y,2 ) ) *sin ( fPhi ) ) /2.;
        n3[2][0]= ( pow ( e1z,2 ) - pow ( e3z,2 ) ) /2. + ( ( pow ( e1z,2 ) + pow ( e3z,2 ) ) *sin ( fPhi ) ) /2.;
        n3[3][0]= ( e1x*e1z - e3x*e3z ) /2. + ( ( e1x*e1z + e3x*e3z ) *sin ( fPhi ) ) /2.;
        n3[4][0]= ( e1y*e1z - e3y*e3z ) /2. + ( ( e1y*e1z + e3y*e3z ) *sin ( fPhi ) ) /2.;
        n3[5][0]= ( e1x*e1y - e3x*e3y ) /2. + ( ( e1x*e1y + e3x*e3y ) *sin ( fPhi ) ) /2.;


        out[0]=n1;
        out[1]=n2;
        out[2]=n3;

        return out;
    }

    NRvector<Doub>  Yields ( NRtensor<Doub> tensor )
    {
        //cout << "tensoes"<<endl;
        //tensor.Print();
        NRvector<Doub> yields ( 3 ),sigma ( 3 );
        NRmatrix<Doub> pt,vec;
        tensor.EigenSystem ( pt,vec );
        sigma[0] = pt[0][0];
        sigma[1] = pt[0][1];
        sigma[2] = pt[0][2];
        //sigma.Print();
        Doub phi1= 1./2.* ( sigma[0]-sigma[2] )+1./2.* ( sigma[0]+sigma[2] ) *sin ( fPhi )-fc*cos ( fPhi );
        Doub phi2= 1./2.* ( sigma[0]-sigma[1] )+1./2.* ( sigma[0]+sigma[1] ) *sin ( fPhi )-fc*cos ( fPhi );
        Doub phi3= 1./2.* ( sigma[1]-sigma[2] )+1./2.* ( sigma[1]+sigma[2] ) *sin ( fPhi )-fc*cos ( fPhi );
        yields[0]=phi1;
        yields[0]=phi2;
        yields[0]=phi3;
        return yields;

    }

    NRvector<Doub>  phi ( NRtensor<Doub> epse )
    {
        NRmatrix<Doub>  tempepsemat, stresstrial,pt,vec;
        NRmatrix<Doub>  C = GetElasticMatrix();
        epse.FromTensorToNRmatrix ( tempepsemat );
        C.Mult ( tempepsemat, stresstrial );
        NRtensor<Doub>  stresstrialtensor;
        epse.FromNRmatrixToTensor ( stresstrial, stresstrialtensor );
        NRvector<Doub> phiv ( 3. );
        //phiv = Yields(stresstrialtensor);

        phiv[0]=PhiInvars ( stresstrialtensor );
        return phiv;
    }

    void updateatributes ( NRvector<MatDoub> mult );

    /**
    * @brief Implements the elastoplastic decomposition whit the Closest Point Projection Method and computes the tangent operator.
     * @param [in] epst total strain tensor
       * @param [in] epsp plastic strain tensor
       * @param [out] projstress  projected stress tensor
       * @param [out] projstrain projected strain tensor
       * @param [out] Dep Elastoplastic tangent operator
       * @param [out] projgamma Plastic multiplier
     */
    void closestpointproj2 ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep, Doub & projgamma );


    /**
    * @brief Implements the elastoplastic decomposition whit the Closest Point Projection Method and computes the tangent operator.
    * @param [in] epst total strain tensor
    * @param [in] epsp plastic strain tensor
    * @param [out] projstress  projected stress tensor
    * @param [out] projstrain projected strain tensor
    * @param [out] Dep Elastoplastic tangent operator
    * @param [out] projgamma Plastic multiplier
    */
    void closestpointproj ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain, NRmatrix<Doub>  & Dep2, Doub & projgamma );

    /**
     * @brief Implements the elastoplastic decomposition whit the Closest Point Projection Method and computes the tangent operator.
     * @param [in] epst total strain tensor
     * @param [in] epsp plastic strain tensor
     * @param [out] projstress  projected stress tensor
     * @param [out] projstrain projected strain tensorr
     */
    bool closestpointproj ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp, NRtensor<Doub>  & projstress, NRtensor<Doub>  & projstrain,Int &whatphi );


    /*void closestpointproj(NRtensor<Doub>  epst, NRtensor<Doub>  epsp,NRtensor<Doub> &projstress, NRtensor<Doub> &projstrain,NRmatrix<Doub> &Dep2, Doub & projgamma)
    {
    	projgamma=0;
    	Dep2.assign(3,3,0.);
    	NRmatrix<Doub>Dep(6,6,0.);
    	srand((unsigned int)time(NULL));
    	float a = 1.e-6;
    	for(Int i=0;i<6;i++)
    	{
    		//if(fabs(epst[i])<1.e-20)epst[i]=1.e-20;
    	}

    	NRtensor<Doub> epstpertub(epst);

    	bool iselastic = closestpointproj(epst, epsp, projstress, projstrain);
    	if(fsetconsistentangent==false)
    	{
    		iselastic=true;
    	}
    	//iselastic=true;
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

    }*/

    void closestpointproj3 ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp,NRtensor<Doub> &projstress, NRtensor<Doub> &projstrain,NRmatrix<Doub> &Dep2, Doub & projgamma )
    {
        projgamma=0;
        Dep2.assign ( 3,3,0. );
        NRmatrix<Doub>Dep ( 6,6,0. );
        Doub a = epst.Norm();
        a*=1.e-6;

        // a=1.e-5;
        NRtensor<Doub> epstpertub ( epst );

        NRtensor<Doub>  epse = epst - epsp;
        NRmatrix<Doub>  C = GetElasticMatrix();

        NRmatrix<Doub>  tempepsemat, stresstrial, Dept;
        epse.FromTensorToNRmatrix ( tempepsemat );

        C.Mult ( tempepsemat, stresstrial );

        NRtensor<Doub>  stresstrialtensor,epsptensor;
        epse.FromNRmatrixToTensor ( stresstrial, stresstrialtensor );

        Int whatphi;
        bool iselastic = closestpointproj ( epst, epsp, projstress, projstrain,whatphi );
        Doub phiinvars=PhiInvars ( stresstrialtensor );
        if ( fsetconsistentangent==false ) {
            //iselastic=true;
        }
        //iselastic=true;
        //iselastic=true;
        if ( iselastic==false ) {

            //2d
            NRvector<Int> Index ( 3 );
            Index[0]=0;//xx
            Index[1]=1;//yy
            Index[2]=5;//xy
            NRvector<NRtensor<Doub>> vecprojperturbstress ( 6 ),vecprojperturbstrain ( 6 );
            for ( Int iperturb=0; iperturb<3; iperturb++ ) {
                Int iit=Index[iperturb];
                NRtensor<Doub> epstpertub ( epst ), projstressperturb, projstrainperturb;
                epstpertub[iit]+=epst[iit]*a;
                //closestpointproj(epst, epsp, projstress, projstrain);
                closestpointproj ( epstpertub, epsp, projstressperturb, projstrainperturb,whatphi );

                NRtensor<Doub> dsig ( projstressperturb ),deps ( epstpertub );
                //NRtensor<Doub> dsig ( projstressperturb ),deps ( projstrainperturb );

                dsig-=projstress;
                deps-=epst;
                //	deps-=projstrain;
                for ( Int ivar=0; ivar<3; ivar++ ) {
                    Int jjt=Index[ivar];
                    Dep2[ivar][iperturb]=dsig[jjt]/deps[iit];
                }

            }
            return;


        } else {
            Dep = GetElasticMatrix();
            Dep2[0][0] = Dep[0][0];
            Dep2[0][1] = Dep[0][1];
            Dep2[0][2] = Dep[0][5];
            Dep2[1][0] = Dep[1][0];
            Dep2[1][1] = Dep[1][1];
            Dep2[1][2] = Dep[1][5];
            Dep2[2][0] = Dep[5][0];
            Dep2[2][1] = Dep[5][1];
            Dep2[2][2] = Dep[5][5];
            return;
        }



    }

    void ComputeNumericalDep2 ( NRtensor<Doub>  epst, NRtensor<Doub>  epsp,NRtensor<Doub> &projstress, NRtensor<Doub> &projstrain,NRmatrix<Doub> &Dep )
    {

        srand ( ( unsigned int ) time ( NULL ) );
        float a = 1.e-8;

        epst.Print();

        NRtensor<Doub> randtensor,projstressperturb,  projstrainperturb;

        randtensor.XX() = ( float ( rand() ) /float ( ( RAND_MAX ) ) );
        randtensor.YY() = ( float ( rand() ) /float ( ( RAND_MAX ) ) );
        randtensor.ZZ() = ( float ( rand() ) /float ( ( RAND_MAX ) ) );
        randtensor.XZ() = ( float ( rand() ) /float ( ( RAND_MAX ) ) );
        randtensor.YZ() = ( float ( rand() ) /float ( ( RAND_MAX ) ) );
        randtensor.XY() = ( float ( rand() ) /float ( ( RAND_MAX ) ) );
        NRtensor<Doub> epstpertub ( epst );
        //for(Int i=0;i<6;i++)epstpertub[i]=epst[i]+randtensor[i] *a;
        //for(Int i=0;i<6;i++)epstpertub[i]=epst[i]+epst[i]*a;
        //epstpertub+=a;
        epstpertub.XY()+=epst.XY() *a;
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

        Int whatphi;
        closestpointproj ( epst, epsp, projstress, projstrain,whatphi );
        cout << "projstress"<<endl;
        projstress.Print();
        closestpointproj ( epstpertub, epsp, projstressperturb, projstrainperturb,whatphi );

        Dep.assign ( 6,6,0. );
        NRtensor<Doub> dsig ( projstressperturb ),deps ( epstpertub );
        cout << "projstressperturb"<<endl;
        projstressperturb.Print();


        dsig-=projstress;
        deps-=epst;

        //for(Int i=0;i<6;i++)deps[i]=a;
        Dep[0][0]=dsig.XX() /deps.XX();
        Dep[0][1]=dsig.XX() /deps.YY();
        Dep[0][2]=dsig.XX() /deps.ZZ();
        Dep[0][3]=dsig.XX() /deps.XZ();
        Dep[0][4]=dsig.XX() /deps.YZ();
        Dep[0][5]=dsig.XX() /deps.XY();

        Dep[1][0]=dsig.YY() /deps.XX();
        Dep[1][1]=dsig.YY() /deps.YY();
        Dep[1][2]=dsig.YY() /deps.ZZ();
        Dep[1][3]=dsig.YY() /deps.XZ();
        Dep[1][4]=dsig.YY() /deps.YZ();
        Dep[1][5]=dsig.YY() /deps.XY();

        Dep[2][0]=dsig.ZZ() /deps.XX();
        Dep[2][1]=dsig.ZZ() /deps.YY();
        Dep[2][2]=dsig.ZZ() /deps.ZZ();
        Dep[2][3]=dsig.ZZ() /deps.XZ();
        Dep[2][4]=dsig.ZZ() /deps.YZ();
        Dep[2][5]=dsig.ZZ() /deps.XY();

        Dep[3][0]=dsig.XZ() /deps.XX();
        Dep[3][1]=dsig.XZ() /deps.YY();
        Dep[3][2]=dsig.XZ() /deps.ZZ();
        Dep[3][3]=dsig.XZ() /deps.XZ();
        Dep[3][4]=dsig.XZ() /deps.YZ();
        Dep[3][5]=dsig.XZ() /deps.XY();

        Dep[4][0]=dsig.YZ() /deps.XX();
        Dep[4][1]=dsig.YZ() /deps.YY();
        Dep[4][2]=dsig.YZ() /deps.ZZ();
        Dep[4][3]=dsig.YZ() /deps.XZ();
        Dep[4][4]=dsig.YZ() /deps.YZ();
        Dep[4][5]=dsig.YZ() /deps.XY();

        Dep[5][0]=dsig.XY() /deps.XX();
        Dep[5][1]=dsig.XY() /deps.YY();
        Dep[5][2]=dsig.XY() /deps.ZZ();
        Dep[5][3]=dsig.XY() /deps.XZ();
        Dep[5][4]=dsig.XY() /deps.YZ();
        Dep[5][5]=dsig.XY() /deps.XY();
    }


    /// Set up the phi
    void SetPhi ( Doub phi )
    {
        fPhi = phi;
    }

    Doub Psi()
    {
        return fPsi;
    }


    Doub Cohesion()
    {
        return fc;
    }

    void SetCohesion ( Doub cohesion )
    {
        fc = cohesion;
    }


    Doub E()
    {
        return fyoung;
    }

    Doub Poisson()
    {
        return fnu;
    }

    void GetMatConstants ( NRvector<Doub>& consts )
    {
        consts.assign ( 5, 0. );
        consts[0] = fyoung;
        consts[1] = fnu;
        consts[2] = fc ;
        consts[3] = fPhi;
        consts[4] = fPsi;

    }

    inline void restoreoriginalatributes()
    {
        SetUp ( fyoung0,fnu0, fc0,fPhi0, fPsi0 );
    }

    void SetMatConstants ( NRvector<Doub>& consts )
    {
		//associativo
        SetUp ( consts[0], consts[1], consts[2], consts[3],consts[3] );
    }

    void SetTangentMatrixType ( bool type )
    {
        fsetconsistentangent = type;
    }
    NRmatrix<Doub>  GetElasticMatrix()
    {
        MatDoub C ( 6, 6, 0. );
        Doub G = fmu, K = fK;
        C[0][0] = ( 4 * G ) / 3 + K;
        C[0][1] = - ( ( 2 * G ) / 3 ) + K;
        C[0][2] = - ( ( 2 * G ) / 3 ) + K;
        C[0][3] = 0.;
        C[0][4] = 0.;
        C[0][5] = 0.;
        C[1][0] = - ( ( 2 * G ) / 3 ) + K;
        C[1][1] = ( 4 * G ) / 3 + K;
        C[1][2] = - ( ( 2 * G ) / 3 ) + K;
        C[1][3] = 0.;
        C[1][4] = 0.;
        C[1][5] = 0.;
        C[2][0] = - ( ( 2 * G ) / 3 ) + K;
        C[2][1] = - ( ( 2 * G ) / 3 ) + K;
        C[2][2] = ( 4 * G ) / 3 + K;
        C[2][3] = 0.;
        C[2][4] = 0.;
        C[2][5] = 0.;
        C[3][0] = 0;
        C[3][1] = 0;
        C[3][2] = 0;
        C[3][3] = G;
        C[3][4] = 0.;
        C[3][5] = 0.;
        C[4][0] = 0;
        C[4][1] = 0;
        C[4][2] = 0;
        C[4][3] = 0.;
        C[4][4] = G;
        C[4][5] = 0.;
        C[5][0] = 0;
        C[5][1] = 0;
        C[5][2] = 0;
        C[5][3] = 0.;
        C[5][4] = 0.;
        C[5][5] = G;
        return C;
    }
    NRmatrix<Doub>   GetInverseElasticMatrix()
    {
        MatDoub C ( 6, 6, 0. );
        Doub G = fmu, K = fK;
        C[0][0] = ( G + 3 * K ) / ( 9.*G*K );
        C[0][1] = -1 / ( 6.*G ) + 1 / ( 9.*K );
        C[0][2] = -1 / ( 6.*G ) + 1 / ( 9.*K );
        C[0][3] = 0.;
        C[0][4] = 0.;
        C[0][5] = 0.;
        C[1][0] = -1 / ( 6.*G ) + 1 / ( 9.*K );
        C[1][1] = ( G + 3 * K ) / ( 9.*G*K );
        C[1][2] = -1 / ( 6.*G ) + 1 / ( 9.*K );
        C[1][3] = 0.;
        C[1][4] = 0.;
        C[1][5] = 0.;
        C[2][0] = -1 / ( 6.*G ) + 1 / ( 9.*K );
        C[2][1] = -1 / ( 6.*G ) + 1 / ( 9.*K );
        C[2][2] = ( G + 3 * K ) / ( 9.*G*K );
        C[2][3] = 0.;
        C[2][4] = 0.;
        C[2][5] = 0.;
        C[3][0] = 0;
        C[3][1] = 0;
        C[3][2] = 0;
        C[3][3] = 1. / G;
        C[3][4] = 0.;
        C[3][5] = 0.;
        C[4][0] = 0;
        C[4][1] = 0;
        C[4][2] = 0;
        C[4][3] = 0.;
        C[4][4] = 1. / G;
        C[4][5] = 0.;
        C[5][0] = 0;
        C[5][1] = 0;
        C[5][2] = 0;
        C[5][3] = 0.;
        C[5][4] = 0.;
        C[5][5] = 1. / G;
        return C;
    }





};


