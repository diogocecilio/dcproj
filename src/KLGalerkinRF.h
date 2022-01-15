#pragma once
#include "elastmat2D.h"
#include "shapequad.h"

#include "eigen_sym.h"

#include "eigen_unsym.h"
#include <random>
#include <cmath>

#include "nr3.h"
#include "elastmat2D.h"
#include "cholesky.h"
#include "mesh.h"
#include "shapequad.h"

class KLGalerkinRF : public elastmat2D
{
public:
    KLGalerkinRF ( Int order, Doub Lx, Doub Ly, Int type, Int samples, Int expansionorder );
    ~KLGalerkinRF();

    void ContributeB ( MatDoub &BE, Doub xi, Doub eta, Doub w, MatDoub elcoords );
    void CacStiffB ( MatDoub &BE, const MatDoub  &elcoords );
    void AssembleB ( MatDoub &B );

    void ContributeC ( MatDoub &CE, MatDoub psis1, MatDoub GradPsi1, MatDoub elcoords1, Doub w1, MatDoub psis2, MatDoub GradPsi2, MatDoub elcoords2, Doub w2 );
    void CacStiffC ( MatDoub &CE, const MatDoub  &elcoords1, const MatDoub  &elcoords2 );
    void AssembleC ( MatDoub &C );

    void SolveGenEigValProblem ( VecComplex & val, MatDoub & vec, NRmatrix<MatDoub> & HHAT, std::vector<std::vector<double>> &errpost );

    void GenerateGaussinRandomField ( VecComplex& val, MatDoub& vec, NRmatrix<MatDoub>& HHAT, std::vector<std::vector<double>>& errpost );

    void GenerateNonGaussinRandomField ( VecComplex& val, MatDoub& vec, NRmatrix<MatDoub>& HHAT, std::vector<std::vector<double>>& errpost );

    Doub AutocorrelationFunc ( MatDoub  x1, MatDoub  x2 );

    Doub PerfomIntegralOfListconst ( const MatDoub &Vec );

    Doub PerfomIntegralOfListconst2 ( const MatDoub &Vec );


    void ComputeVarianceError ( VecComplex &val, MatDoub &vec, MatDoub &error );

    void OutPutPost ( std::vector<std::vector<double>> & postdata, std::ofstream &file )
    {
        file.clear();
        for ( Int i = 0; i < postdata.size(); i++ ) {
            for ( Int j = 0; j < postdata[0].size(); j++ ) {
                file << postdata[i][j] << " ";
            }
            file << endl;
        }
        file.close();
    }

    inline Doub GetLx()
    {
        return fLx;
    }
    inline Doub GetLy()
    {
        return fLy;
    }
    inline Doub GetSig()
    {
        return 0.;
    }
    inline Doub Getftype()
    {
        return ftype;
    }
    inline Doub GetExpansionorder()
    {
        return fexpansionorder;
    }
    inline Doub GetSamples()
    {
        return fsamples;
    }

    inline void SetMesh ( mesh *inmesh )
    {
        fmesh = *inmesh;
    }

private:
    Doub fyoung;
    Doub fnu;
    Doub fbodyforce;
    Int fplanestress;
    Doub fthickness;
    Int fOrder;
    Doub fLx;
    Doub fLy;
    //Doub fsig;
    Int ftype;
    Int fsamples;
    Int fexpansionorder;
    mesh fmesh;
	shape *fshape;
};

