//
// Created by Diogo Cecílio on 10/12/21.
//

#pragma once
#include "nr3.h"
class shape
{
public:

    shape();
    ~shape();

    virtual void shapes ( MatDoub&psis, MatDoub &gradpsis, Doub xi, Doub eta ) =0;
    virtual void shapes ( MatrixXd &psis, MatrixXd &gradpsis, double xi, double eta ) =0;
    virtual void shapes1D ( MatDoub&psis, MatDoub &gradpsis, Doub xi ) =0;
    virtual void pointsandweigths ( MatDoub & pts ) =0;
    virtual void pointsandweigths1D ( MatDoub & pts ) =0;
    virtual NRmatrix<Doub> GetBaseNodes() =0;
};


