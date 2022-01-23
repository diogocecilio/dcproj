#pragma once
#include "nr3.h"
#include "shape.h"

class shapetri: public shape
{
public:
    shapetri ( int order, int type, Int dim=2 )
    {
        forder = order;
        ftype = type;
    }

    shapetri()
    {

    }
    ~shapetri()
    {

    }

    void shapes ( MatDoub&psis, MatDoub &gradpsis, Doub xi, Doub eta )
    {
        switch ( forder ) {
        case 1:

            psis.resize ( 3, 1 );

            gradpsis.resize ( 2, 3 );

            psis[0][0] = 1 - eta - xi;
            psis[1][0] = xi;
            psis[2][0] = eta;

            gradpsis[0][0] = -1;
            gradpsis[0][1] = 1.;
            gradpsis[0][2] = 0.;

            gradpsis[1][0] =-1.;
            gradpsis[1][1] = 0.;
            gradpsis[1][2] = 1.;
			
			
// 			psis[0][0] = eta;
//             psis[1][0] = xi;
//             psis[2][0] = 1 - eta - xi;
// 
//             gradpsis[0][0] = 0;
//             gradpsis[0][1] = 1.;
//             gradpsis[0][2] = -1.;
// 
//             gradpsis[1][0] =1.;
//             gradpsis[1][1] = 0.;
//             gradpsis[1][2] = -1.;
			
// 			psis[0][0] = xi;
//             psis[1][0] = eta;
//             psis[2][0] = 1 - eta - xi;
// 
//             gradpsis[0][0] = 1.;
//             gradpsis[0][1] = 0.;
//             gradpsis[0][2] = -1.;
// 
//             gradpsis[1][0] =0.;
//             gradpsis[1][1] = 1.;
//             gradpsis[1][2] = -1.;
			
/*			
			psis[0][0] = xi;
            psis[1][0] = 1 - eta - xi;
            psis[2][0] = eta;

            gradpsis[0][0] = 1.;
            gradpsis[0][1] = -1.;
            gradpsis[0][2] = 0.;

            gradpsis[1][0] =0.;
            gradpsis[1][1] = -1.;
            gradpsis[1][2] = 1.;*/
			

            break;

        case 2:
            psis.assign ( 6, 1, 0.0000 );
            gradpsis.assign ( 2, 6, 0.000 );


			
			psis[0][0] = (-1. + eta + xi)*(-1. + 2.*eta + 2.*xi);
            psis[1][0] =xi*(-1. + 2.*xi);
            psis[2][0] = eta*(-1. + 2.*eta);
            psis[3][0] = -4.*xi*(-1. + eta + xi);
            psis[4][0] = 4.*eta*xi;
            psis[5][0] =-4.*eta*(-1. + eta + xi);


            gradpsis[0][0] =-3. + 4.*eta + 4.*xi;
            gradpsis[0][1] = -1. + 4.*xi;
            gradpsis[0][2] = 0.;
            gradpsis[0][3] = -4.*(-1. + eta + 2.*xi);
            gradpsis[0][4] = 4.*eta;
            gradpsis[0][5] = -4.*eta;

			
			gradpsis[1][0] = -3. + 4.*eta + 4.*xi;
            gradpsis[1][1] = 0.;
            gradpsis[1][2] =-1. + 4.*eta;
            gradpsis[1][3] =-4.*xi;
            gradpsis[1][4] = 4.*xi;
            gradpsis[1][5] = -4.*(-1. + 2.*eta + xi);
			
// 			psis[2][0] = (-1. + eta + xi)*(-1. + 2.*eta + 2.*xi);
//             psis[1][0] =xi*(-1. + 2.*xi);
//             psis[0][0] = eta*(-1. + 2.*eta);
//             psis[5][0] = -4.*xi*(-1. + eta + xi);
//             psis[4][0] = 4.*eta*xi;
//             psis[3][0] =-4.*eta*(-1. + eta + xi);
// 
// 
//             gradpsis[0][2] =-3. + 4.*eta + 4.*xi;
//             gradpsis[0][1] = -1. + 4.*xi;
//             gradpsis[0][0] = 0.;
//             gradpsis[0][5] = -4.*(-1. + eta + 2.*xi);
//             gradpsis[0][4] = 4.*eta;
//             gradpsis[0][3] = -4.*eta;
// 
// 			
// 			gradpsis[1][2] = -3. + 4.*eta + 4.*xi;
//             gradpsis[1][1] = 0.;
//             gradpsis[1][0] =-1. + 4.*eta;
//             gradpsis[1][5] =-4.*xi;
//             gradpsis[1][4] = 4.*xi;
//             gradpsis[1][3] = -4.*(-1. + 2.*eta + xi);
			


            break;
        }

    }


    void shapes ( MatrixXd &psis, MatrixXd &gradpsis, double xi, double eta )
    {
        DebugStop();
    }



    void shapes1D ( MatDoub&psis, MatDoub &gradpsis, Doub xi )
    {
        switch ( forder ) {
        case 1:
            psis.resize ( 2, 1 );
            gradpsis.resize ( 2, 1 );
            psis[0][0] = 0.5 - 0.5*xi;
            psis[1][0] = 0.5 * ( 1. + xi );

            gradpsis[0][0] = ( 1 - xi ) / 2;
            gradpsis[1][0] = ( 1 + xi ) / 2;
            break;
        case 2:
            psis.resize ( 3, 1 );
            psis[0][0] = -0.5 * xi + 0.5 * xi*xi;
            psis[1][0] = 1. - 1. * xi *xi;
            psis[2][0] = 0.5 * xi * ( 1. + xi );


            //List(((-1 + xi)*xi) / 2., 1 - Power(xi, 2), (xi*(1 + xi)) / 2.)

            gradpsis.resize ( 3, 1 );
            gradpsis[0][0] = -0.5 + xi;
            gradpsis[1][0] = -2.* xi;
            gradpsis[2][0] = 0.5 + xi;
            break;
        case 3:
            //{-0.0625 + 0.0625 xi + 0.5625 xi ^ 2 - 0.5625 xi ^ 3,
            //	0.5625 - 1.6875 xi - 0.5625 xi ^ 2 +
            //	1.6875 xi ^ 3, -1.6875 (0.333333 + xi) (1. + xi) (-1. + 1. xi),
            //	0.5625 (-0.333333 + xi) (0.333333 + xi) (1. + xi)}
            break;
        }

    }

    void pointsandweigths ( MatDoub & pts )
    {

        if ( true ) {
            if ( forder == 1 ) {


                double co[1][3] = {
                    { 1./3.,1./3.,1./2. }
                };
                pts.resize ( 1, 3 );
                for ( Int i = 0; i < 1; i++ ) {
                    for ( Int j = 0; j < 3; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }
                

            } else {

//                 double co[3][3] = {{1./6., 1./6., 1./2. *1./3.}, {2./3., 1./6., 1./2. *1./3.}, {1./6., 2./3., 1./2. *1./3.}};
// 
//                 pts.resize ( 3, 3 );
//                 for ( Int i = 0; i < 3; i++ ) {
//                     for ( Int j = 0; j < 3; j++ ) {
//                         pts[i][j] = co[i][j];
//                     }
//                 }
                


//                 double co[4][3] ={{1./3., 1./3., -27./48.*1./2.}, {1./5., 1./5., 25./48.* 1./2.}, {1./5., 3./5., 
//  25./48.* 1./2.}, {3./5., 1./5., 25./48.* 1./2.}};
// 
//                 pts.resize ( 4, 3 );
//                 for ( Int i = 0; i < 4; i++ ) {
//                     for ( Int j = 0; j < 3; j++ ) {
//                         pts[i][j] = co[i][j];
//                     }
//                 }

Doub r1=0.1012865073235;
Doub r2=0.7974269853531;
Doub r3=r1;
Doub r4=0.4701420641051;
Doub r5=r4;
Doub r6=0.0597158717898;
Doub r7=1./3.;
Doub s1=r1;
Doub s2=r1;
Doub s3=r2;
Doub s4=r6;
Doub s5=r4;
Doub s6=r4;
Doub s7 = r7;
Doub w1=0.1259391805448;
Doub w2=w1;
Doub w3=w1;
Doub w4=0.1323941527885;
Doub w5=w4;
Doub w6=w4;
Doub w7=0.225;
                 double co[7][3] ={
					 {r1,s1, w1}, 
					 {r2,s2, w2},
					 {r3,s3, w3},
					 {r4, s4, w4},
					 {r5, s5, w5},
					 {r6, s6, w6},
					 {r7, s7, w7}};
                pts.resize ( 7, 3 );
                for ( Int i = 0; i < 7; i++ ) {
                    for ( Int j = 0; j < 3; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }
            }


        }

    }

    void pointsandweigths1D ( MatDoub & pts )
    {
        //double co[4][2] = { { -0.861136, 0.347855 },{ -0.339981, 0.652145 },{ 
//0.339981,
        //		0.652145 },{ 0.861136, 0.347855 } };

        double co[7][2] = { {-0.949108, 0.129485}, { -0.741531, 0.279705 }, {
                -0.405845,
                    0.38183
                }, { 0., 0.417959 }, { 0.405845, 0.38183 }, {
                0.741531,
                0.279705
            }, { 0.949108, 0.129485 }
        };

        pts.resize ( 7, 2 );
        for ( Int i = 0; i < pts.nrows(); i++ ) {
            for ( Int j = 0; j < 2; j++ ) {
                pts[i][j] = co[i][j];
            }
        }
    }


    NRmatrix<Doub> GetBaseNodes()
    {
        NRmatrix<Doub> pts;
        if ( forder == 2 ) {
            pts.assign ( 6,2,0. );

            //double v[6][2] = { {1.,0.},{0.,1.},{0.,0.},{1/2.,1./2.},
            //     {0.,1./2.},{1./2.,0.}
            // };
            double v[6][2] = { {0.,0.},{1.,0.},{0.,1.},{1./2.,0.},
                {1./2.,1./2.},{0.,1./2.}
            };
            for ( Int i =0; i<6; i++ ) {
                pts[i][0]=v[i][0];
                pts[i][1]=v[i][1];
            }

        } else if ( forder == 1 ) {
            pts.assign ( 3,2,0. );
            //(0,1)
            //------
            //---------
            //-----------
            //------------
            //-------------
            //(0,0)-----------(1,0)
            double v[3][2] = {
                {0,0},{1,0},{0,1}
            };
            for ( Int i =0; i<3; i++ ) {
                pts[i][0]=v[i][0];
                pts[i][1]=v[i][1];
            }
        }
        return pts;
    }

protected:
    int forder;
    int ftype;

};


