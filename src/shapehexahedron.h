#pragma once
#include "nr3.h"
#include "shape.h"

class shapehexahedron
{
public:

    shapehexahedron ( int order, int type,Int dim=3 )
    {
        forder = order;
        ftype = type;
    }

    shapehexahedron()
    {
    }
    ~shapehexahedron()
    {

    }

    void shapes ( MatDoub&psis, MatDoub &gradpsis, Doub xi, Doub eta, Doub zeta )
    {
        switch ( forder ) {
        case 1:
            psis.assign ( 8, 1, 0. );
            gradpsis.assign ( 3, 8, 0. );

            //psis
            psis[0][0] = ( ( 1 - eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8.;
            psis[1][0] = ( ( 1 - eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8.;
            psis[2][0] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8.;
            psis[3][0] = ( ( 1 + eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8.;
            psis[4][0] = ( ( 1 - eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8.;
            psis[5][0] = ( ( 1 - eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8.;
            psis[6][0] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8.;
            psis[7][0] = ( ( 1 + eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8.;

            //dpsi/dxi
            gradpsis[0][0] = - ( ( 1 - eta ) * ( 1 - zeta ) ) /8.;
            gradpsis[0][1] = ( ( 1 - eta ) * ( 1 - zeta ) ) /8.;
            gradpsis[0][2] = ( ( 1 + eta ) * ( 1 - zeta ) ) /8.;
            gradpsis[0][3] = - ( ( 1 + eta ) * ( 1 - zeta ) ) /8.;
            gradpsis[0][4] = - ( ( 1 - eta ) * ( 1 + zeta ) ) /8.;
            gradpsis[0][5] = ( ( 1 - eta ) * ( 1 + zeta ) ) /8.;
            gradpsis[0][6] = ( ( 1 + eta ) * ( 1 + zeta ) ) /8.;
            gradpsis[0][7] = - ( ( 1 + eta ) * ( 1 + zeta ) ) /8.;

            //dpsi/deta
            gradpsis[1][0] = - ( ( 1 - xi ) * ( 1 - zeta ) ) /8.;
            gradpsis[1][1] = - ( ( 1 + xi ) * ( 1 - zeta ) ) /8.;
            gradpsis[1][2] = ( ( 1 + xi ) * ( 1 - zeta ) ) /8.;
            gradpsis[1][3] = ( ( 1 - xi ) * ( 1 - zeta ) ) /8.;
            gradpsis[1][4] = - ( ( 1 - xi ) * ( 1 + zeta ) ) /8.;
            gradpsis[1][5] = - ( ( 1 + xi ) * ( 1 + zeta ) ) /8.;
            gradpsis[1][6] = ( ( 1 + xi ) * ( 1 + zeta ) ) /8.;
            gradpsis[1][7] = ( ( 1 - xi ) * ( 1 + zeta ) ) /8.;

            //dpsi/zeta

            gradpsis[2][0] = - ( ( 1 - eta ) * ( 1 - xi ) ) /8.;
            gradpsis[2][1] = - ( ( 1 - eta ) * ( 1 + xi ) ) /8.;
            gradpsis[2][2] = - ( ( 1 + eta ) * ( 1 + xi ) ) /8.;
            gradpsis[2][3] = - ( ( 1 + eta ) * ( 1 - xi ) ) /8.;
            gradpsis[2][4] = ( ( 1 - eta ) * ( 1 - xi ) ) /8.;
            gradpsis[2][5] = ( ( 1 - eta ) * ( 1 + xi ) ) /8.;
            gradpsis[2][6] = ( ( 1 + eta ) * ( 1 + xi ) ) /8.;
            gradpsis[2][7] = ( ( 1 + eta ) * ( 1 - xi ) ) /8.;


            break;
        case 2:
            psis.assign ( 20, 1, 0. );
            gradpsis.assign ( 3, 20, 0. );

            //psis
            psis[0][0] = - ( ( 1 - eta ) * ( 1 - xi ) * ( 1 - zeta ) * ( 2 + eta + xi + zeta ) ) /8.;
            psis[1][0] = - ( ( 1 - eta ) * ( 1 + xi ) * ( 1 - zeta ) * ( 2 + eta - xi + zeta ) ) /8.;
            psis[2][0] = - ( ( 1 + eta ) * ( 1 + xi ) * ( 1 - zeta ) * ( 2 - eta - xi + zeta ) ) /8.;
            psis[3][0] = - ( ( 1 + eta ) * ( 1 - xi ) * ( 1 - zeta ) * ( 2 - eta + xi + zeta ) ) /8.;

            psis[4][0] = - ( ( 1 - eta ) * ( 1 - xi ) * ( 2 + eta + xi - zeta ) * ( 1 + zeta ) ) /8.;
            psis[5][0] = - ( ( 1 - eta ) * ( 1 + xi ) * ( 2 + eta - xi - zeta ) * ( 1 + zeta ) ) /8.;
            psis[6][0] = - ( ( 1 + eta ) * ( 1 + xi ) * ( 2 - eta - xi - zeta ) * ( 1 + zeta ) ) /8.;
            psis[7][0] = - ( ( 1 + eta ) * ( 1 - xi ) * ( 2 - eta + xi - zeta ) * ( 1 + zeta ) ) /8.;

            psis[8][0] = ( ( 1 - eta ) * ( 1 - xi*xi ) * ( 1 - zeta ) ) /4.;
            psis[9][0] = ( ( 1 - eta*eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /4.;
            psis[10][0] = ( ( 1 + eta ) * ( 1 - xi*xi ) * ( 1 - zeta ) ) /4.;
            psis[11][0] = ( ( 1 - eta*eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /4.;


            //GID numbering
            psis[16][0] = ( ( 1 - eta ) * ( 1 - xi*xi ) * ( 1 + zeta ) ) /4.;
            psis[17][0] = ( ( 1 - eta*eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /4.;
            psis[18][0] = ( ( 1 + eta ) * ( 1 - xi*xi ) * ( 1 + zeta ) ) /4.;
            psis[19][0] = ( ( 1 - eta*eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /4.;


            psis[12][0] = ( ( 1 - eta ) * ( 1 - xi ) * ( 1 - zeta*zeta ) ) /4.;
            psis[13][0] = ( ( 1 - eta ) * ( 1 + xi ) * ( 1 - zeta*zeta ) ) /4.;
            psis[14][0] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 - zeta*zeta ) ) /4.;
            psis[15][0] = ( ( 1 + eta ) * ( 1 - xi ) * ( 1 - zeta*zeta ) ) /4.;


            //dpsi/dxi
            gradpsis[0][0] = - ( ( 1 - eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8. + ( ( 1 - eta ) * ( 1 - zeta ) * ( 2 + eta + xi + zeta ) ) /8.;
            gradpsis[0][1] = ( ( 1 - eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8. - ( ( 1 - eta ) * ( 1 - zeta ) * ( 2 + eta - xi + zeta ) ) /8.;
            gradpsis[0][2] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8. - ( ( 1 + eta ) * ( 1 - zeta ) * ( 2 - eta - xi + zeta ) ) /8.;
            gradpsis[0][3] = - ( ( 1 + eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8. + ( ( 1 + eta ) * ( 1 - zeta ) * ( 2 - eta + xi + zeta ) ) /8.;

            gradpsis[0][4] = - ( ( 1 - eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8. + ( ( 1 - eta ) * ( 2 + eta + xi - zeta ) * ( 1 + zeta ) ) /8.;
            gradpsis[0][5] = ( ( 1 - eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8. - ( ( 1 - eta ) * ( 2 + eta - xi - zeta ) * ( 1 + zeta ) ) /8.;
            gradpsis[0][6] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8. - ( ( 1 + eta ) * ( 2 - eta - xi - zeta ) * ( 1 + zeta ) ) /8.;
            gradpsis[0][7] = - ( ( 1 + eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8. + ( ( 1 + eta ) * ( 2 - eta + xi - zeta ) * ( 1 + zeta ) ) /8.;

            gradpsis[0][8] = - ( ( 1 - eta ) *xi* ( 1 - zeta ) ) /2.;
            gradpsis[0][9] = ( ( 1 - eta*eta ) * ( 1 - zeta ) ) /4.;
            gradpsis[0][10] = - ( ( 1 + eta ) *xi* ( 1 - zeta ) ) /2.;
            gradpsis[0][11] = - ( ( 1 - eta*eta ) * ( 1 - zeta ) ) /4.;

            gradpsis[0][16] = - ( ( 1 - eta ) *xi* ( 1 + zeta ) ) /2.;
            gradpsis[0][17] = ( ( 1 - eta*eta ) * ( 1 + zeta ) ) /4.;
            gradpsis[0][18] = - ( ( 1 + eta ) *xi* ( 1 + zeta ) ) /2.;
            gradpsis[0][19] = - ( ( 1 - eta*eta ) * ( 1 + zeta ) ) /4.;

            gradpsis[0][12] = - ( ( 1 - eta ) * ( 1 - zeta*zeta ) ) /4.;
            gradpsis[0][13] = ( ( 1 - eta ) * ( 1 - zeta*zeta ) ) /4.;
            gradpsis[0][14] = ( ( 1 + eta ) * ( 1 - zeta*zeta ) ) /4.;
            gradpsis[0][15] = - ( ( 1 + eta ) * ( 1 - zeta*zeta ) ) /4.;

            //dpsi/deta

            gradpsis[1][0] =  - ( ( 1 - eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8. + ( ( 1 - xi ) * ( 1 - zeta ) * ( 2 + eta + xi + zeta ) ) /8.;
            gradpsis[1][1] =  - ( ( 1 - eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8. + ( ( 1 + xi ) * ( 1 - zeta ) * ( 2 + eta - xi + zeta ) ) /8.;
            gradpsis[1][2] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8. - ( ( 1 + xi ) * ( 1 - zeta ) * ( 2 - eta - xi + zeta ) ) /8.;
            gradpsis[1][3] = ( ( 1 + eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8. - ( ( 1 - xi ) * ( 1 - zeta ) * ( 2 - eta + xi + zeta ) ) /8.;

            gradpsis[1][4] =  - ( ( 1 - eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8. + ( ( 1 - xi ) * ( 2 + eta + xi - zeta ) * ( 1 + zeta ) ) /8.;
            gradpsis[1][5] =  - ( ( 1 - eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8. + ( ( 1 + xi ) * ( 2 + eta - xi - zeta ) * ( 1 + zeta ) ) /8.;
            gradpsis[1][6] = ( ( 1 + eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8. - ( ( 1 + xi ) * ( 2 - eta - xi - zeta ) * ( 1 + zeta ) ) /8.;
            gradpsis[1][7] = ( ( 1 + eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8. - ( ( 1 - xi ) * ( 2 - eta + xi - zeta ) * ( 1 + zeta ) ) /8.;

            gradpsis[1][8] =  - ( ( 1 - xi*xi ) * ( 1 - zeta ) ) /4.;
            gradpsis[1][9] =  - ( eta* ( 1 + xi ) * ( 1 - zeta ) ) /2.;
            gradpsis[1][10] = ( ( 1 - xi*xi ) * ( 1 - zeta ) ) /4.;
            gradpsis[1][11] = - ( eta* ( 1 - xi ) * ( 1 - zeta ) ) /2.;

            gradpsis[1][16] = - ( ( 1 - xi*xi ) * ( 1 + zeta ) ) /4.;
            gradpsis[1][17] = - ( eta* ( 1 + xi ) * ( 1 + zeta ) ) /2.;
            gradpsis[1][18] = ( ( 1 - xi*xi ) * ( 1 + zeta ) ) /4.;
            gradpsis[1][19] = - ( eta* ( 1 - xi ) * ( 1 + zeta ) ) /2.;

            gradpsis[1][12] = - ( ( 1 - xi ) * ( 1 - zeta*zeta ) ) /4.;
            gradpsis[1][13] = - ( ( 1 + xi ) * ( 1 - zeta*zeta ) ) /4.;
            gradpsis[1][14] = ( ( 1 + xi ) * ( 1 - zeta*zeta ) ) /4.;
            gradpsis[1][15] = ( ( 1 - xi ) * ( 1 - zeta*zeta ) ) /4.;

            //dpsi/zeta

            gradpsis[2][0] =  - ( ( 1 - eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8. + ( ( 1 - eta ) * ( 1 - xi ) * ( 2 + eta + xi + zeta ) ) /8.;
            gradpsis[2][1] =  - ( ( 1 - eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8. + ( ( 1 - eta ) * ( 1 + xi ) * ( 2 + eta - xi + zeta ) ) /8.;
            gradpsis[2][2] =  - ( ( 1 + eta ) * ( 1 + xi ) * ( 1 - zeta ) ) /8. + ( ( 1 + eta ) * ( 1 + xi ) * ( 2 - eta - xi + zeta ) ) /8.;
            gradpsis[2][3] =  - ( ( 1 + eta ) * ( 1 - xi ) * ( 1 - zeta ) ) /8. + ( ( 1 + eta ) * ( 1 - xi ) * ( 2 - eta + xi + zeta ) ) /8.;

            gradpsis[2][4] =  - ( ( 1 - eta ) * ( 1 - xi ) * ( 2 + eta + xi - zeta ) ) /8. + ( ( 1 - eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8.;
            gradpsis[2][5] =  - ( ( 1 - eta ) * ( 1 + xi ) * ( 2 + eta - xi - zeta ) ) /8. + ( ( 1 - eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8.;
            gradpsis[2][6] =  - ( ( 1 + eta ) * ( 1 + xi ) * ( 2 - eta - xi - zeta ) ) /8. + ( ( 1 + eta ) * ( 1 + xi ) * ( 1 + zeta ) ) /8.;
            gradpsis[2][7] =  - ( ( 1 + eta ) * ( 1 - xi ) * ( 2 - eta + xi - zeta ) ) /8. + ( ( 1 + eta ) * ( 1 - xi ) * ( 1 + zeta ) ) /8.;

            gradpsis[2][8] =  - ( ( 1 - eta ) * ( 1 - xi*xi ) ) /4.;
            gradpsis[2][9] =  - ( ( 1 - eta*eta ) * ( 1 + xi ) ) /4.;
            gradpsis[2][10] = - ( ( 1 + eta ) * ( 1 - xi*xi ) ) /4.;
            gradpsis[2][11] = - ( ( 1 - eta*eta ) * ( 1 - xi ) ) /4.;

            gradpsis[2][16] = ( ( 1 - eta ) * ( 1 - xi*xi ) ) /4.;
            gradpsis[2][17] = ( ( 1 - eta*eta ) * ( 1 + xi ) ) /4.;
            gradpsis[2][18] = ( ( 1 + eta ) * ( 1 - xi*xi ) ) /4.;
            gradpsis[2][19] = ( ( 1 - eta*eta ) * ( 1 - xi ) ) /4.;

            gradpsis[2][12] =  - ( ( 1 - eta ) * ( 1 - xi ) *zeta ) /2.;
            gradpsis[2][13] =  - ( ( 1 - eta ) * ( 1 + xi ) *zeta ) /2.;
            gradpsis[2][14] =  - ( ( 1 + eta ) * ( 1 + xi ) *zeta ) /2.;
            gradpsis[2][15] =  - ( ( 1 + eta ) * ( 1 - xi ) *zeta ) /2.;

            break;
        }

    }

    /*void shapes(MatDoub&psis, MatDoub &gradpsis, Doub xi, Doub eta, Doub zeta)
    {
    	switch (forder)
    	{
    	case 1:
            psis.assign(8, 1, 0.);
    		gradpsis.assign(3, 8, 0.);

            //psis
    		psis[0][0] = ((1 - eta)*(1 - xi)*(1 - zeta))/8.;
    		psis[1][0] = ((1 - eta)*(1 + xi)*(1 - zeta))/8.;
    		psis[2][0] = ((1 + eta)*(1 + xi)*(1 - zeta))/8.;
    		psis[3][0] = ((1 + eta)*(1 - xi)*(1 - zeta))/8.;
            psis[4][0] = ((1 - eta)*(1 - xi)*(1 + zeta))/8.;
    		psis[5][0] = ((1 - eta)*(1 + xi)*(1 + zeta))/8.;
    		psis[6][0] = ((1 + eta)*(1 + xi)*(1 + zeta))/8.;
    		psis[7][0] = ((1 + eta)*(1 - xi)*(1 + zeta))/8.;

            //dpsi/dxi
            gradpsis[0][0] = -((1 - eta)*(1 - zeta))/8.;
            gradpsis[0][1] =  ((1 - eta)*(1 - zeta))/8.;
            gradpsis[0][2] =  ((1 + eta)*(1 - zeta))/8.;
            gradpsis[0][3] = -((1 + eta)*(1 - zeta))/8.;
            gradpsis[0][4] = -((1 - eta)*(1 + zeta))/8.;
            gradpsis[0][5] =  ((1 - eta)*(1 + zeta))/8.;
            gradpsis[0][6] =  ((1 + eta)*(1 + zeta))/8.;
            gradpsis[0][7] = -((1 + eta)*(1 + zeta))/8.;

            //dpsi/deta
            gradpsis[1][0] = -((1 - xi)*(1 - zeta))/8.;
            gradpsis[1][1] = -((1 + xi)*(1 - zeta))/8.;
            gradpsis[1][2] =  ((1 + xi)*(1 - zeta))/8.;
            gradpsis[1][3] =  ((1 - xi)*(1 - zeta))/8.;
            gradpsis[1][4] = -((1 - xi)*(1 + zeta))/8.;
            gradpsis[1][5] = -((1 + xi)*(1 + zeta))/8.;
            gradpsis[1][6] =  ((1 + xi)*(1 + zeta))/8.;
            gradpsis[1][7] =  ((1 - xi)*(1 + zeta))/8.;

            //dpsi/zeta

            gradpsis[2][0] = -((1 - eta)*(1 - xi))/8.;
            gradpsis[2][1] = -((1 - eta)*(1 + xi))/8.;
            gradpsis[2][2] = -((1 + eta)*(1 + xi))/8.;
            gradpsis[2][3] = -((1 + eta)*(1 - xi))/8.;
            gradpsis[2][4] =  ((1 - eta)*(1 - xi))/8.;
            gradpsis[2][5] =  ((1 - eta)*(1 + xi))/8.;
            gradpsis[2][6] =  ((1 + eta)*(1 + xi))/8.;
            gradpsis[2][7] =  ((1 + eta)*(1 - xi))/8.;


    		break;
    	case 2:
    		psis.assign(20, 1, 0.);
    		gradpsis.assign(3, 20, 0.);

            //psis
            psis[0][0] = -((1 - eta)*(1 - xi)*(1 - zeta)*(2 + eta + xi + zeta))/8.;
            psis[1][0] = -((1 - eta)*(1 + xi)*(1 - zeta)*(2 + eta - xi + zeta))/8.;
            psis[2][0] = -((1 + eta)*(1 + xi)*(1 - zeta)*(2 - eta - xi + zeta))/8.;
            psis[3][0] = -((1 + eta)*(1 - xi)*(1 - zeta)*(2 - eta + xi + zeta))/8.;

            psis[4][0] = -((1 - eta)*(1 - xi)*(2 + eta + xi - zeta)*(1 + zeta))/8.;
            psis[5][0] = -((1 - eta)*(1 + xi)*(2 + eta - xi - zeta)*(1 + zeta))/8.;
            psis[6][0] = -((1 + eta)*(1 + xi)*(2 - eta - xi - zeta)*(1 + zeta))/8.;
            psis[7][0] = -((1 + eta)*(1 - xi)*(2 - eta + xi - zeta)*(1 + zeta))/8.;

            psis[8][0] =  ((1 - eta)*(1 - xi*xi)*(1 - zeta))/4.;
            psis[9][0] =  ((1 - eta*eta)*(1 + xi)*(1 - zeta))/4.;
            psis[10][0] = ((1 + eta)*(1 - xi*xi)*(1 - zeta))/4.;
            psis[11][0] = ((1 - eta*eta)*(1 - xi)*(1 - zeta))/4.;

            psis[12][0] = ((1 - eta)*(1 - xi*xi)*(1 + zeta))/4.;
            psis[13][0] = ((1 - eta*eta)*(1 + xi)*(1 + zeta))/4.;
            psis[14][0] = ((1 + eta)*(1 - xi*xi)*(1 + zeta))/4.;
            psis[15][0] = ((1 - eta*eta)*(1 - xi)*(1 + zeta))/4.;


            psis[16][0] = ((1 - eta)*(1 - xi)*(1 - zeta*zeta))/4.;
            psis[17][0] = ((1 - eta)*(1 + xi)*(1 - zeta*zeta))/4.;
            psis[18][0] = ((1 + eta)*(1 + xi)*(1 - zeta*zeta))/4.;
            psis[19][0] = ((1 + eta)*(1 - xi)*(1 - zeta*zeta))/4.;


            //dpsi/dxi
            gradpsis[0][0] = -((1 - eta)*(1 - xi)*(1 - zeta))/8. + ((1 - eta)*(1 - zeta)*(2 + eta + xi + zeta))/8.;
            gradpsis[0][1] =  ((1 - eta)*(1 + xi)*(1 - zeta))/8. - ((1 - eta)*(1 - zeta)*(2 + eta - xi + zeta))/8.;
            gradpsis[0][2] =  ((1 + eta)*(1 + xi)*(1 - zeta))/8. - ((1 + eta)*(1 - zeta)*(2 - eta - xi + zeta))/8.;
            gradpsis[0][3] = -((1 + eta)*(1 - xi)*(1 - zeta))/8. + ((1 + eta)*(1 - zeta)*(2 - eta + xi + zeta))/8.;

            gradpsis[0][4] = -((1 - eta)*(1 - xi)*(1 + zeta))/8. + ((1 - eta)*(2 + eta + xi - zeta)*(1 + zeta))/8.;
            gradpsis[0][5] =  ((1 - eta)*(1 + xi)*(1 + zeta))/8. - ((1 - eta)*(2 + eta - xi - zeta)*(1 + zeta))/8.;
            gradpsis[0][6] =  ((1 + eta)*(1 + xi)*(1 + zeta))/8. - ((1 + eta)*(2 - eta - xi - zeta)*(1 + zeta))/8.;
            gradpsis[0][7] = -((1 + eta)*(1 - xi)*(1 + zeta))/8. + ((1 + eta)*(2 - eta + xi - zeta)*(1 + zeta))/8.;

            gradpsis[0][8] = -((1 - eta)*xi*(1 - zeta))/2.;
            gradpsis[0][9] =  ((1 - eta*eta)*(1 - zeta))/4.;
            gradpsis[0][10] = -((1 + eta)*xi*(1 - zeta))/2.;
            gradpsis[0][11] = -((1 - eta*eta)*(1 - zeta))/4.;

            gradpsis[0][12] = -((1 - eta)*xi*(1 + zeta))/2.;
            gradpsis[0][13] =  ((1 - eta*eta)*(1 + zeta))/4.;
            gradpsis[0][14] = -((1 + eta)*xi*(1 + zeta))/2.;
            gradpsis[0][15] = -((1 - eta*eta)*(1 + zeta))/4.;

            gradpsis[0][16] = -((1 - eta)*(1 - zeta*zeta))/4.;
            gradpsis[0][17] =  ((1 - eta)*(1 - zeta*zeta))/4.;
            gradpsis[0][18] =  ((1 + eta)*(1 - zeta*zeta))/4.;
            gradpsis[0][19] = -((1 + eta)*(1 - zeta*zeta))/4.;

            //dpsi/deta

            gradpsis[1][0] =  -((1 - eta)*(1 - xi)*(1 - zeta))/8. + ((1 - xi)*(1 - zeta)*(2 + eta + xi + zeta))/8.;
            gradpsis[1][1] =  -((1 - eta)*(1 + xi)*(1 - zeta))/8. + ((1 + xi)*(1 - zeta)*(2 + eta - xi + zeta))/8.;
            gradpsis[1][2] =   ((1 + eta)*(1 + xi)*(1 - zeta))/8. - ((1 + xi)*(1 - zeta)*(2 - eta - xi + zeta))/8.;
            gradpsis[1][3] =   ((1 + eta)*(1 - xi)*(1 - zeta))/8. - ((1 - xi)*(1 - zeta)*(2 - eta + xi + zeta))/8.;

            gradpsis[1][4] =  -((1 - eta)*(1 - xi)*(1 + zeta))/8. + ((1 - xi)*(2 + eta + xi - zeta)*(1 + zeta))/8.;
            gradpsis[1][5] =  -((1 - eta)*(1 + xi)*(1 + zeta))/8. + ((1 + xi)*(2 + eta - xi - zeta)*(1 + zeta))/8.;
            gradpsis[1][6] =   ((1 + eta)*(1 + xi)*(1 + zeta))/8. - ((1 + xi)*(2 - eta - xi - zeta)*(1 + zeta))/8.;
            gradpsis[1][7] =   ((1 + eta)*(1 - xi)*(1 + zeta))/8. - ((1 - xi)*(2 - eta + xi - zeta)*(1 + zeta))/8.;

            gradpsis[1][8] =  -((1 - xi*xi)*(1 - zeta))/4.;
            gradpsis[1][9] =  -(eta*(1 + xi)*(1 - zeta))/2.;
            gradpsis[1][10] =  ((1 - xi*xi)*(1 - zeta))/4.;
            gradpsis[1][11] = -(eta*(1 - xi)*(1 - zeta))/2.;

            gradpsis[1][12] = -((1 - xi*xi)*(1 + zeta))/4.;
            gradpsis[1][13] = -(eta*(1 + xi)*(1 + zeta))/2.;
            gradpsis[1][14] =  ((1 - xi*xi)*(1 + zeta))/4.;
            gradpsis[1][15] = -(eta*(1 - xi)*(1 + zeta))/2.;

            gradpsis[1][16] = -((1 - xi)*(1 - zeta*zeta))/4.;
            gradpsis[1][17] = -((1 + xi)*(1 - zeta*zeta))/4.;
            gradpsis[1][18] =  ((1 + xi)*(1 - zeta*zeta))/4.;
            gradpsis[1][19] =  ((1 - xi)*(1 - zeta*zeta))/4.;

            //dpsi/zeta

            gradpsis[2][0] =  -((1 - eta)*(1 - xi)*(1 - zeta))/8. + ((1 - eta)*(1 - xi)*(2 + eta + xi + zeta))/8.;
            gradpsis[2][1] =  -((1 - eta)*(1 + xi)*(1 - zeta))/8. + ((1 - eta)*(1 + xi)*(2 + eta - xi + zeta))/8.;
            gradpsis[2][2] =  -((1 + eta)*(1 + xi)*(1 - zeta))/8. + ((1 + eta)*(1 + xi)*(2 - eta - xi + zeta))/8.;
            gradpsis[2][3] =  -((1 + eta)*(1 - xi)*(1 - zeta))/8. + ((1 + eta)*(1 - xi)*(2 - eta + xi + zeta))/8.;

            gradpsis[2][4] =  -((1 - eta)*(1 - xi)*(2 + eta + xi - zeta))/8. + ((1 - eta)*(1 - xi)*(1 + zeta))/8.;
            gradpsis[2][5] =  -((1 - eta)*(1 + xi)*(2 + eta - xi - zeta))/8. + ((1 - eta)*(1 + xi)*(1 + zeta))/8.;
            gradpsis[2][6] =  -((1 + eta)*(1 + xi)*(2 - eta - xi - zeta))/8. + ((1 + eta)*(1 + xi)*(1 + zeta))/8.;
            gradpsis[2][7] =  -((1 + eta)*(1 - xi)*(2 - eta + xi - zeta))/8. + ((1 + eta)*(1 - xi)*(1 + zeta))/8.;

            gradpsis[2][8] =  -((1 - eta)*(1 - xi*xi))/4.;
            gradpsis[2][9] =  -((1 - eta*eta)*(1 + xi))/4.;
            gradpsis[2][10] = -((1 + eta)*(1 - xi*xi))/4.;
            gradpsis[2][11] = -((1 - eta*eta)*(1 - xi))/4.;

            gradpsis[2][12] =  ((1 - eta)*(1 - xi*xi))/4.;
            gradpsis[2][13] =  ((1 - eta*eta)*(1 + xi))/4.;
            gradpsis[2][14] =  ((1 + eta)*(1 - xi*xi))/4.;
            gradpsis[2][15] =  ((1 - eta*eta)*(1 - xi))/4.;

            gradpsis[2][16] =  -((1 - eta)*(1 - xi)*zeta)/2.;
            gradpsis[2][17] =  -((1 - eta)*(1 + xi)*zeta)/2.;
            gradpsis[2][18] =  -((1 + eta)*(1 + xi)*zeta)/2.;
            gradpsis[2][19] =  -((1 + eta)*(1 - xi)*zeta)/2.;


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
                double co[9][4] = {
                    {   1., 0., 0., 1.3333333333},
                    {-1., 0., 0., 1.3333333333},
                    {  0., 1., 0., 1.3333333333},
                    { 0., -1., 0., 1.3333333333},
                    {  0., 0., 1., 1.3333333333},
                    { 0., 0., -1., 1.3333333333}
                };
                pts.resize ( 6, 4 );
                for ( Int i = 0; i < 6; i++ ) {
                    for ( Int j = 0; j < 4; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }
            } else {

                double co[14][4] = {
                    {0.7958224257, 0.0000000000, 0.0000000000, 0.8864265927},
                    {-0.7958224257, 0.0000000000, 0.0000000000, 0.8864265927},
                    {0.0000000000, 0.7958224257, 0.0000000000, 0.8864265927},
                    {0.0000000000, -0.7958224257, 0.0000000000, 0.8864265927},
                    {0.0000000000, 0.0000000000, 0.7958224257, 0.8864265927},
                    {0.0000000000, 0.0000000000, -0.7958224257, 0.8864265927},
                    {0.7587869106, 0.7587869106, 0.7587869106, 0.3351800554},
                    {0.7587869106, -0.7587869106, 0.7587869106, 0.3351800554},
                    {0.7587869106, 0.7587869106, -0.7587869106, 0.3351800554},
                    {0.7587869106, -0.7587869106, -0.7587869106, 0.3351800554},
                    {-0.7587869106, 0.7587869106, 0.7587869106, 0.3351800554},
                    {-0.7587869106, -0.7587869106, 0.7587869106, 0.3351800554},
                    {-0.7587869106, 0.7587869106, -0.7587869106, 0.3351800554},
                    {-0.7587869106, -0.7587869106, -0.7587869106, 0.3351800554}
                };
                pts.resize ( 14, 4 );
                for ( Int i = 0; i < 14; i++ ) {
                    for ( Int j = 0; j < 4; j++ ) {
                        pts[i][j] = co[i][j];
                    }
                }
            }


        }







        //double co[36][3] = { {-0.93247, 0.93247, 0.0293521}, {-0.661209, 0.93247,
        //0.0618073}, {-0.238619, 0.93247, 0.0801651}, {0.238619, 0.93247,
        //0.0801651}, {0.661209, 0.93247, 0.0618073}, {0.93247, 0.93247,
        //0.0293521}, {-0.93247, 0.661209, 0.0618073}, {-0.661209, 0.661209,
        //0.130149}, {-0.238619, 0.661209, 0.168805}, {0.238619, 0.661209,
        //0.168805}, {0.661209, 0.661209, 0.130149}, {0.93247, 0.661209,
        //0.0618073}, {-0.93247, 0.238619, 0.0801651}, {-0.661209, 0.238619,
        //0.168805}, {-0.238619, 0.238619, 0.218943}, {0.238619, 0.238619,
        //0.218943}, {0.661209, 0.238619, 0.168805}, {0.93247, 0.238619,
        //0.0801651}, {-0.93247, -0.238619, 0.0801651}, {-0.661209, -0.238619,
        // 0.168805}, {-0.238619, -0.238619, 0.218943}, {0.238619, -0.238619,
        //0.218943}, {0.661209, -0.238619, 0.168805}, {0.93247, -0.238619,
        //0.0801651}, {-0.93247, -0.661209, 0.0618073}, {-0.661209, -0.661209,
        // 0.130149}, {-0.238619, -0.661209, 0.168805}, {0.238619, -0.661209,
        //0.168805}, {0.661209, -0.661209, 0.130149}, {0.93247, -0.661209,
        //0.0618073}, {-0.93247, -0.93247, 0.0293521}, {-0.661209, -0.93247,
        //0.0618073}, {-0.238619, -0.93247, 0.0801651}, {0.238619, -0.93247,
        //0.0801651}, {0.661209, -0.93247, 0.0618073}, {0.93247, -0.93247,
        //0.0293521} };

        //pts.resize(36, 3);
        //for (Int i = 0; i < 36; i++)
        //{
        //	for (Int j = 0; j < 3; j++)
        //	{
        //		pts[i][j] = co[i][j];
        //	}
        //}


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



protected:
    int forder;
    int ftype;

};


