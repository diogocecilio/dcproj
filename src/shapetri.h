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
			
            break;
			
        case 2:
            psis.assign ( 6, 1, 0.0000 );
            gradpsis.assign ( 2, 6, 0.000 );

            psis[0][0] = 2.*(0.5 - 1.*eta - 1.*xi)*(1. - 1.*eta - 1.*xi);
            psis[1][0] =2.*(-0.5 + xi)*xi;
            psis[2][0] = 2.*(-0.5 + eta)*eta;
            psis[3][0] = 4.*(1. - 1.*eta - 1.*xi)*xi;
            psis[4][0] = 4.*eta*xi;
            psis[5][0] =4.*eta*(1. - 1.*eta - 1.*xi);


            gradpsis[0][0] =-2.*(0.5 - 1.*eta - 1.*xi) - 2.*(1. - 1.*eta - 1.*xi);
            gradpsis[0][1] = 2.*(-0.5 + xi) + 2.*xi;
            gradpsis[0][2] = 0.;
            gradpsis[0][3] = 4.*(1. - 1.*eta - 1.*xi) - 4.*xi;
            gradpsis[0][4] = 4.*eta;
            gradpsis[0][5] = -4.*eta;
 

            gradpsis[1][0] = -2.*(0.5 - 1.*eta - 1.*xi) - 2.*(1. - 1.*eta - 1.*xi);
            gradpsis[1][1] = 0.;
            gradpsis[1][2] =2.*(-0.5 + eta) + 2.*eta;
            gradpsis[1][3] =-4.*xi;
            gradpsis[1][4] = 4.*xi;
            gradpsis[1][5] = -4.*eta + 4.*(1. - 1.*eta - 1.*xi);

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
                    { 1./3.,1./3.,0.5 }
                };
                pts.resize ( 1, 3 );
                for ( Int i = 0; i < 1; i++ ) {
                    for ( Int j = 0; j < 3; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }

            } else {

                double co[3][3] = {
				{0.166667, 0.166667, 0.166667}, 
				{0.666667, 0.166667,0.166667}, 
				{0.166667, 0.666667, 0.166667}
				};
                
                pts.resize ( 3, 3 );
                for ( Int i = 0; i < 3; i++ ) {
                    for ( Int j = 0; j < 3; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }

            }


        }

    }

    void pointsandweigths1D ( MatDoub & pts )
    {
        //double co[4][2] = { { -0.861136, 0.347855 },{ -0.339981, 0.652145 },{ 0.339981,
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
		pts.assign( 3,2,0. );
            //(0,1)
            //------
            //---------
            //-----------
            //------------
            //-------------
            //(0,0)-----------(1,0)
            double v[3][2] = {
                {1,0},{0,1},{0,0}
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


