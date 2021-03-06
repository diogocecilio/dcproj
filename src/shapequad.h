#pragma once
#include "nr3.h"
#include "shape.h"

class shapequad: public shape
{
public:
    shapequad ( int order, int type, Int dim=2 )
    {
        forder = order;
        ftype = type;
    }

    shapequad()
    {

    }
    ~shapequad()
    {

    }

    void shapes ( MatDoub&psis, MatDoub &gradpsis, Doub xi, Doub eta )
    {
        switch ( forder ) {
        case 1:
            psis.resize ( 4, 1 );
            gradpsis.resize ( 2, 4 );
            psis[0][0] = ( ( 1 - eta ) * ( 1 - xi ) ) / 4.;
            psis[1][0] = ( ( 1 - eta ) * ( 1 + xi ) ) / 4.;
            psis[2][0] = ( ( 1 + eta ) * ( 1 + xi ) ) / 4.;
            psis[3][0] = ( ( 1 + eta ) * ( 1 - xi ) ) / 4.;

            gradpsis[0][0] = 1. / 4. * ( -1. + eta );
            gradpsis[0][1] = ( 1. - eta ) / 4.;
            gradpsis[0][2] = ( 1. + eta ) / 4.;
            gradpsis[0][3] = 1. / 4. * ( -1. - eta );

            gradpsis[1][0] = 1. / 4. * ( -1. + xi );
            gradpsis[1][1] = 1. / 4. * ( -1. - xi );
            gradpsis[1][2] = ( 1. + xi ) / 4.;
            gradpsis[1][3] = ( 1. - xi ) / 4.;
            break;
        case 2:
            psis.assign ( 8, 1, 0.0000 );
            gradpsis.assign ( 2, 8, 0.000 );

            psis[0][0] = 1. / 4. * ( 1. - eta ) * ( 1. - xi ) * ( -1. - eta - xi );
            psis[1][0] = 1. / 4. * ( 1. - eta ) * ( 1. + xi ) * ( -1. - eta + xi );
            psis[2][0] = 1. / 4.* ( 1. + eta ) * ( 1. + xi ) * ( -1. + eta + xi );
            psis[3][0] = 1. / 4. * ( 1. + eta ) * ( 1. - xi ) * ( -1. + eta - xi );
            psis[4][0] = 1. / 2. * ( 1. - eta ) * ( 1 - xi *xi );
            psis[5][0] = 1. / 2. * ( 1. - eta *eta ) * ( 1. + xi );
            psis[6][0] = 1. / 2. * ( 1. + eta ) * ( 1. - xi *xi );
            psis[7][0] = 1. / 2. * ( 1. - eta *eta ) * ( 1. - xi );


            gradpsis[0][0] = - ( ( 1 - eta ) * ( 1 - xi ) ) / 4. - ( ( 1 - eta ) * ( -1 - eta - xi ) ) / 4.;
            gradpsis[0][1] = ( ( 1 - eta ) * ( 1 + xi ) ) / 4. + ( ( 1 - eta ) * ( -1 - eta + xi ) ) / 4.;
            gradpsis[0][2] = ( ( 1 + eta ) * ( 1 + xi ) ) / 4. + ( ( 1 + eta ) * ( -1 + eta + xi ) ) / 4.;
            gradpsis[0][3] = - ( ( 1 + eta ) * ( 1 - xi ) ) / 4. - ( ( 1 + eta ) * ( -1 + eta - xi ) ) / 4.;
            gradpsis[0][4] = ( -1 + eta ) *xi;
            gradpsis[0][5] = ( 1 - eta*eta ) / 2.;
            gradpsis[0][6] = ( -1 - eta ) *xi;
            gradpsis[0][7] = ( -1 + pow ( eta, 2 ) ) / 2.;

            gradpsis[1][0] = - ( ( 1 - eta ) * ( 1 - xi ) ) / 4. - ( ( 1 - xi ) * ( -1 - eta - xi ) ) / 4.;
            gradpsis[1][1] = - ( ( 1 - eta ) * ( 1 + xi ) ) / 4. - ( ( 1 + xi ) * ( -1 - eta + xi ) ) / 4.;
            gradpsis[1][2] = ( ( 1 + eta ) * ( 1 + xi ) ) / 4. + ( ( 1 + xi ) * ( -1 + eta + xi ) ) / 4.;
            gradpsis[1][3] = ( ( 1 + eta ) * ( 1 - xi ) ) / 4. + ( ( 1 - xi ) * ( -1 + eta - xi ) ) / 4.;
            gradpsis[1][4] = ( -1 + xi*xi ) / 2.;
            gradpsis[1][5] = -eta* ( 1. + xi );
            gradpsis[1][6] = 1. / 2. * ( 1. - xi *xi );
            gradpsis[1][7] = -eta* ( 1. - xi );
            break;
        }

    }
 /*   
        void shapes( MatDoub&psis, MatDoub &gradpsis, Doub xi, Doub eta )
    {
        switch ( forder ) {
        case 1:
            psis.resize ( 4, 1 );
            gradpsis.resize ( 2, 4 );
            psis[0][0] = ( ( 1 - eta ) * ( 1 - xi ) ) / 4.;
            psis[1][0] = ( ( 1 - eta ) * ( 1 + xi ) ) / 4.;
            psis[2][0] = ( ( 1 + eta ) * ( 1 + xi ) ) / 4.;
            psis[3][0] = ( ( 1 + eta ) * ( 1 - xi ) ) / 4.;

            gradpsis[0][0] = 1. / 4. * ( -1. + eta );
            gradpsis[0][1] = ( 1. - eta ) / 4.;
            gradpsis[0][2] = ( 1. + eta ) / 4.;
            gradpsis[0][3] = 1. / 4. * ( -1. - eta );

            gradpsis[1][0] = 1. / 4. * ( -1. + xi );
            gradpsis[1][1] = 1. / 4. * ( -1. - xi );
            gradpsis[1][2] = ( 1. + xi ) / 4.;
            gradpsis[1][3] = ( 1. - xi ) / 4.;
            break;
        case 2:
            psis.assign ( 8, 1, 0.0000 );
            gradpsis.assign ( 2, 8, 0.000 );

            psis[0][0] =-((-1 + eta)*(-1 + xi)*(1 + eta + xi))/4.;
            psis[1][0] =((-1 + eta)*(-1 + pow(xi,2)))/2.;
            psis[2][0] =((-1 + eta)*(1 + eta - xi)*(1 + xi))/4.;
            psis[3][0] =-((-1 + pow(eta,2))*(1 + xi))/2.;
            psis[4][0] =  ((1 + eta)*(1 + xi)*(-1 + eta + xi))/4.;
            psis[5][0] = -((1 + eta)*(-1 + pow(xi,2)))/2.;
            psis[6][0] = ((-1 + xi)*(1 - pow(eta,2) + xi + eta*xi))/4.;
            psis[7][0] =((-1 + pow(eta,2))*(-1 + xi))/2.;


            gradpsis[0][0] =-((-1 + eta)*(-1 + xi))/4. - ((-1 + eta)*(1 + eta + xi))/4.;
            gradpsis[0][1] =(-1 + eta)*xi ;
            gradpsis[0][2] = ((-1 + eta)*(1 + eta - xi))/4. - ((-1 + eta)*(1 + xi))/4.;
            gradpsis[0][3] = (1 - pow(eta,2))/2.;
            gradpsis[0][4] =  ((1 + eta)*(1 + xi))/4. + ((1 + eta)*(-1 + eta + xi))/4.;
            gradpsis[0][5] = -((1 + eta)*xi);
            gradpsis[0][6] = ((1 + eta)*(-1 + xi))/4. + (1 - pow(eta,2) + xi + eta*xi)/4.;
            gradpsis[0][7] = (-1 + pow(eta,2))/2.;

            gradpsis[1][0] = -((-1 + eta)*(-1 + xi))/4. - ((-1 + xi)*(1 + eta + xi))/4.;
            gradpsis[1][1] =  (-1 + pow(xi,2))/2.;
            gradpsis[1][2] =  ((-1 + eta)*(1 + xi))/4. + ((1 + eta - xi)*(1 + xi))/4.;
            gradpsis[1][3] =  -(eta*(1 + xi));
            gradpsis[1][4] =  ((1 + eta)*(1 + xi))/4. + ((1 + xi)*(-1 + eta + xi))/4.;
            gradpsis[1][5] =  (1 - pow(xi,2))/2.;
            gradpsis[1][6] =  ((-1 + xi)*(-2*eta + xi))/4.;
            gradpsis[1][7] = eta*(-1 + xi);
            break;
        }

    }*/


    void shapes ( MatrixXd &psis, MatrixXd &gradpsis, double xi, double eta )
    {
        switch ( forder ) {
        case 1:
            psis.resize ( 4, 1 );
            gradpsis.resize ( 2, 4 );
            psis ( 0,0 ) = ( ( 1 - eta ) * ( 1 - xi ) ) / 4.;
            psis ( 1,0 ) = ( ( 1 - eta ) * ( 1 + xi ) ) / 4.;
            psis ( 2,0 ) = ( ( 1 + eta ) * ( 1 + xi ) ) / 4.;
            psis ( 3,0 ) = ( ( 1 + eta ) * ( 1 - xi ) ) / 4.;

            gradpsis ( 0,0 ) = 1. / 4. * ( -1. + eta );
            gradpsis ( 0,1 ) = ( 1. - eta ) / 4.;
            gradpsis ( 0,2 ) = ( 1. + eta ) / 4.;
            gradpsis ( 0,3 ) = 1. / 4. * ( -1. - eta );

            gradpsis ( 1,0 ) = 1. / 4. * ( -1. + xi );
            gradpsis ( 1,1 ) = 1. / 4. * ( -1. - xi );
            gradpsis ( 1,2 ) = ( 1. + xi ) / 4.;
            gradpsis ( 1,3 ) = ( 1. - xi ) / 4.;
            break;
        case 2:
            psis.resize ( 8, 1 );
            gradpsis.resize ( 2, 8 );

            psis ( 0,0 ) = 1. / 4. * ( 1. - eta ) * ( 1. - xi ) * ( -1. - eta - xi );
            psis ( 1,0 ) = 1. / 4. * ( 1. - eta ) * ( 1. + xi ) * ( -1. - eta + xi );
            psis ( 2,0 ) = 1. / 4.* ( 1. + eta ) * ( 1. + xi ) * ( -1. + eta + xi );
            psis ( 3,0 ) = 1. / 4. * ( 1. + eta ) * ( 1. - xi ) * ( -1. + eta - xi );
            psis ( 4,0 ) = 1. / 2. * ( 1. - eta ) * ( 1 - xi *xi );
            psis ( 5,0 ) = 1. / 2. * ( 1. - eta *eta ) * ( 1. + xi );
            psis ( 6,0 ) = 1. / 2. * ( 1. + eta ) * ( 1. - xi *xi );
            psis ( 7,0 ) = 1. / 2. * ( 1. - eta *eta ) * ( 1. - xi );


            gradpsis ( 0,0 ) = - ( ( 1 - eta ) * ( 1 - xi ) ) / 4. - ( ( 1 - eta ) * ( -1 - eta - xi ) ) / 4.;
            gradpsis ( 0,1 ) = ( ( 1 - eta ) * ( 1 + xi ) ) / 4. + ( ( 1 - eta ) * ( -1 - eta + xi ) ) / 4.;
            gradpsis ( 0,2 ) = ( ( 1 + eta ) * ( 1 + xi ) ) / 4. + ( ( 1 + eta ) * ( -1 + eta + xi ) ) / 4.;
            gradpsis ( 0,3 ) = - ( ( 1 + eta ) * ( 1 - xi ) ) / 4. - ( ( 1 + eta ) * ( -1 + eta - xi ) ) / 4.;
            gradpsis ( 0,4 ) = ( -1 + eta ) *xi;
            gradpsis ( 0,5 ) = ( 1 - eta*eta ) / 2.;
            gradpsis ( 0,6 ) = ( -1 - eta ) *xi;
            gradpsis ( 0,7 ) = ( -1 + pow ( eta, 2 ) ) / 2.;

            gradpsis ( 1,0 ) = - ( ( 1 - eta ) * ( 1 - xi ) ) / 4. - ( ( 1 - xi ) * ( -1 - eta - xi ) ) / 4.;
            gradpsis ( 1,1 ) = - ( ( 1 - eta ) * ( 1 + xi ) ) / 4. - ( ( 1 + xi ) * ( -1 - eta + xi ) ) / 4.;
            gradpsis ( 1,2 ) = ( ( 1 + eta ) * ( 1 + xi ) ) / 4. + ( ( 1 + xi ) * ( -1 + eta + xi ) ) / 4.;
            gradpsis ( 1,3 ) = ( ( 1 + eta ) * ( 1 - xi ) ) / 4. + ( ( 1 - xi ) * ( -1 + eta - xi ) ) / 4.;
            gradpsis ( 1,4 ) = ( -1 + xi*xi ) / 2.;
            gradpsis ( 1,5 ) = -eta* ( 1. + xi );
            gradpsis ( 1,6 ) = 1. / 2. * ( 1. - xi *xi );
            gradpsis ( 1,7 ) = -eta* ( 1. - xi );
            break;
        }

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


                double co[9][3] = { { -0.774597, 0.774597, 0.308642 },{ 0., 0.774597, 0.493827 },{
                        0.774597,
                        0.774597, 0.308642
                    },{ -0.774597, 0., 0.493827 },{
                        0., 0.,
                        0.790123
                    },{ 0.774597, 0., 0.493827 },{
                        -0.774597, -0.774597,
                            0.308642
                        },{ 0., -0.774597, 0.493827 },{
                        0.774597, -0.774597,
                        0.308642
                    }
                };
                pts.resize ( 9, 3 );
                for ( Int i = 0; i < 9; i++ ) {
                    for ( Int j = 0; j < 3; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }
            } else {

                double co[9][3] = {
                    { -0.774597, 0.774597, 0.308642 },{ 0., 0.774597, 0.493827 },{
                        0.774597,
                        0.774597, 0.308642
                    },{ -0.774597, 0., 0.493827 },{
                        0., 0.,
                        0.790123
                    },{ 0.774597, 0., 0.493827 },{
                        -0.774597, -0.774597,
                            0.308642
                        },{ 0., -0.774597, 0.493827 },{
                        0.774597, -0.774597,
                        0.308642
                    }
                };
                pts.resize ( 9, 3 );
                for ( Int i = 0; i < 9; i++ ) {
                    for ( Int j = 0; j < 3; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }
                /*double co[16][3] = { {-0.861136, 0.861136, 0.121003}, { -0.339981, 0.861136,
                	0.226852 }, { 0.339981, 0.861136, 0.226852 }, { 0.861136, 0.861136,
                	0.121003 }, { -0.861136, 0.339981, 0.226852 }, { -0.339981, 0.339981,
                	0.425293 }, { 0.339981, 0.339981, 0.425293 }, { 0.861136, 0.339981,
                	0.226852 }, { -0.861136, -0.339981, 0.226852 }, { -0.339981, -0.339981,
                	0.425293 }, { 0.339981, -0.339981, 0.425293 }, { 0.861136, -0.339981,
                	0.226852 }, { -0.861136, -0.861136, 0.121003 }, { -0.339981, -0.861136,
                	0.226852 }, { 0.339981, -0.861136, 0.226852 }, { 0.861136, -0.861136,
                	0.121003 } };
                pts.resize(16, 3);
                for (Int i = 0; i < 16; i++)
                {
                	for (Int j = 0; j < 3; j++)
                	{
                		pts[i][j] = co[i][j];
                	}
                }*/

                
//                 double co[36][3] = { {-0.93247, 0.93247, 0.0293521}, {-0.661209, 0.93247,
//                 0.0618073}, {-0.238619, 0.93247, 0.0801651}, {0.238619, 0.93247,
//                 0.0801651}, {0.661209, 0.93247, 0.0618073}, {0.93247, 0.93247,
//                 0.0293521}, {-0.93247, 0.661209, 0.0618073}, {-0.661209, 0.661209,
//                 0.130149}, {-0.238619, 0.661209, 0.168805}, {0.238619, 0.661209,
//                 0.168805}, {0.661209, 0.661209, 0.130149}, {0.93247, 0.661209,
//                 0.0618073}, {-0.93247, 0.238619, 0.0801651}, {-0.661209, 0.238619,
//                 0.168805}, {-0.238619, 0.238619, 0.218943}, {0.238619, 0.238619,
//                 0.218943}, {0.661209, 0.238619, 0.168805}, {0.93247, 0.238619,
//                 0.0801651}, {-0.93247, -0.238619, 0.0801651}, {-0.661209, -0.238619,
//                 0.168805}, {-0.238619, -0.238619, 0.218943}, {0.238619, -0.238619,
//                 0.218943}, {0.661209, -0.238619, 0.168805}, {0.93247, -0.238619,
//                 0.0801651}, {-0.93247, -0.661209, 0.0618073}, {-0.661209, -0.661209,
//                 0.130149}, {-0.238619, -0.661209, 0.168805}, {0.238619, -0.661209,
//                 0.168805}, {0.661209, -0.661209, 0.130149}, {0.93247, -0.661209,
//                 0.0618073}, {-0.93247, -0.93247, 0.0293521}, {-0.661209, -0.93247,
//                 0.0618073}, {-0.238619, -0.93247, 0.0801651}, {0.238619, -0.93247,
//                 0.0801651}, {0.661209, -0.93247, 0.0618073}, {0.93247, -0.93247,
//                 0.0293521} };

//                 pts.resize(36, 3);
//                 for (Int i = 0; i < 36; i++)
//                 {
//                     for (Int j = 0; j < 3; j++)
//                     {
//                         pts[i][j] = co[i][j];
//                     }
//                 }

                         

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
        NRmatrix<Doub> pts ( 8,2,0. );
        if ( forder == 2 ) {
            //(2)----(5)--(1)
            //-------------
            //-------------
            //(6)---------(4)
            //-------------
            //-------------
            //(3)---(7)---(0)
//             double v[8][2] ={{1, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}};

            double v[8][2] = { {-1.,-1.},{1.,-1.},{1.,1.},{-1.,1.},
                {0.,-1.},{1.,0.},{0.,1.},{-1.,0.}
            };
            for ( Int i =0; i<8; i++ ) {
                pts[i][0]=v[i][0];
                pts[i][1]=v[i][1];
            }

        } else if ( forder == 1 ) {

            //(-1,1)-------------(1,1)
            //-----------------------
            //-----------------------
            //------------------------
            //-----------------------
            //-----------------------
            //(-1,-1)-----------(1,-1)
            double v[4][2] = {
                {-1,-1},{1,-1},{1,1},{-1,1},
            };
            for ( Int i =0; i<4; i++ ) {
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


