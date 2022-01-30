#pragma once
#pragma once


#include "nr3.h"


//#include "roots.h"

#include "mesh.h"
#include <math.h>
#include <cmath>

#include "elastoplasticbase.h"
//#include<Eigen/SparseCholesky>
#include <iostream>
//#include <Eigen/Dense>

//#include <direct.h>
#include "elastoplastic2D.h"
#include "KLGalerkinRF.h"
#include "vtkmesh.h"
//using Eigen::MatrixXd;
using namespace Eigen;
class slopeproject
{
public:
    slopeproject ( mesh* inmesh,	KLGalerkinRF * inklgalerking );
    slopeproject ( mesh* inmesh, KLGalerkinRF* inklgalerking, NRmatrix<MatDoub> randomfield );
    slopeproject ( mesh inmesh, KLGalerkinRF inklgalerking, NRmatrix<MatDoub> randomfield );
    slopeproject();
    ~slopeproject();
    slopeproject ( slopeproject& copy );
    void CreateRandomField ( string namefolder );
    void PrintMCS ( string namefolder,int imc,bool print );
    void OutPutPost ( MatDoub& postdata, std::ofstream& file );
    void OutPutPost ( MatInt& postdata, std::ofstream& file );
    void OutPutPost ( std::vector<std::vector<double>>& postdata, std::ofstream& file );
    void OutPutFile ( MatDoub& postdata, std::ofstream& file );
    void OutPutFile1var ( MatDoub& postdata, std::ofstream& file );
    void OutPutFile4var ( MatDoub& postdata, std::ofstream& file );
    void PrintMathematicaFormat ( MatDoub postdata, std::ofstream& file );
    Doub computelamda ( MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l );
    Doub computelamda0 ( MatDoub& dwb,  MatDoub& dw, Doub& l);
    void OutPutPost2 ( std::vector<std::vector<double>>& postdata, std::ofstream& file );
    void ReadMesh ( std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord );
    std::vector<Doub>   vecstr_to_vecdoub2 ( std::vector<string> vs );
    template <class T>
    std::vector<T>   vecstr_to_vec ( std::vector<string> vs );
    void  ReadMatDoub ( MatDoub& matdoub, std::string  file );
    std::vector<Doub>   vecstr_to_vecdoub ( std::vector<string> vs );
    std::vector<Int>   vecstr_to_vecint ( std::vector<string> vs );
    std::vector<std::vector<double>> IterativeProcessNew ( Int ndesi, Doub dlamb0, Doub maxlfac, Int imc );
    std::vector<std::vector<double>> IterativeProcess ( int ndesi, Doub dlamb0, Doub alphatol, int niter );
    std::vector<std::vector<double>> IterativeProcessArcLengthSRM ( int ndesi, Doub dlamb0, Doub alphatol, int niter );
    void   IterativeProcess2();
    std::vector<std::vector<double>>    IterativeProcessGIMBinarySearch();
    void   IterativeProcessSRMBinarySearch();
    void SolveEigenSparse ( int type,MatDoub A, MatDoub b, MatDoub& x );
    void SolveEigen ( SparseMatrix<double> A, VectorXd b, VectorXd& x );
    void SolveEigen ( MatDoub A, MatDoub b, MatDoub& x );
    void SolveEigen2 ( MatDoub A, MatDoub b, MatDoub& x );
    void SolveEigen3 ( MatDoub A, MatDoub b, MatDoub& x );
    void InserBC ( MatDoub& KG, MatDoub& R, MatDoub& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat );
    void InserBC ( SparseMatrix<double> & KG, VectorXd& R, VectorXd& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material*mat );
    void GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub& elcoords );
    void FindIdsInPath ( const MatDoub& path, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idpath );
    void Line ( VecDoub a, VecDoub b, Int ndivs, MatDoub& path );

    void findbcids ( mesh* gmesh, std::vector<std::vector<int>>& idsvector );
	void findbcidsfat ( mesh* gmesh, std::vector<std::vector<int>>& idsvector );
	void findbcidsfat5022 ( mesh* gmesh, std::vector<std::vector<int>>& idsvector );
	void findbcids45( mesh* gmesh, std::vector<std::vector<int>>& idsvector );
	void findbcidsfatsn ( mesh* gmesh, std::vector<std::vector<int>>& idsvector );
    MatDoub  AssembleHhationho ( Int i );
    std::vector<std::vector<double>>  IterativeProcessShearRed ( Doub fac, Doub delta,Doub tol );
    std::vector<std::vector<double>>  IterativeProcessShearRed2 ( Doub fac, Doub delta,Doub tol );
    void MonteCarloSRM ( int iter, int iter2, bool print, string writenamefolder );
    void MonteCarloGIM ( int iter, int iter2,bool print, string writenamefolder );
    mesh* fmesh;
    KLGalerkinRF* fklgalerking;
    NRmatrix<MatDoub> frandomfield;
	std::vector<std::vector<int>> fidsvector;
    void PostVtk (  Int step  );

    template <class T>
    Doub rtsafe ( T &funcd, const Doub x1, const Doub x2, const Doub xacc )
    {
        const Int MAXIT = 1000;
        Doub xh, xl;
        Doub fl = funcd ( x1 );
        Doub fh = funcd ( x2 );
        if ( ( fl > 0.0 && fh > 0.0 ) || ( fl < 0.0 && fh < 0.0 ) )
            std::cout<< "Root must be bracketed in rtsafe" << std::endl;
        if ( fl == 0.0 ) return x1;
        if ( fh == 0.0 ) return x2;
        if ( fl < 0.0 ) {
            xl = x1;
            xh = x2;
        } else {
            xh = x1;
            xl = x2;
        }
        Doub rts = 0.5* ( x1 + x2 );
        Doub dxold = abs ( x2 - x1 );
        Doub dx = dxold;
        Doub f = funcd ( rts );
        Doub df = funcd.df ( rts );
        for ( Int j = 0; j<MAXIT; j++ ) {
            if ( ( ( ( rts - xh ) *df - f ) * ( ( rts - xl ) *df - f ) > 0.0 )
                    || ( abs ( 2.0*f ) > abs ( dxold*df ) ) ) {
                dxold = dx;
                dx = 0.5* ( xh - xl );
                rts = xl + dx;
                if ( xl == rts ) return rts;
            } else {
                dxold = dx;
                dx = f / df;
                Doub temp = rts;
                rts -= dx;
                if ( temp == rts ) return rts;
            }
            if ( abs ( dx ) < xacc ) return rts;
            Doub f = funcd ( rts );
            Doub df = funcd.df ( rts );
            if ( f < 0.0 )
                xl = rts;
            else
                xh = rts;
        }

    }


    NRvector<Doub> computelamda2 ( MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l )
    {
        Int sz = dwb.nrows();
        Doub aa = 0.;
        for ( Int i = 0; i < sz; i++ ) aa += dwb[i][0] * dwb[i][0];
        Doub bb = 0.;
        MatDoub dwcopy = dw;
        dwcopy += dws;
        for ( Int i = 0; i < sz; i++ ) bb += dwb[i][0] * dwcopy[i][0];
        bb *= 2;
        Doub cc = 0.;
        for ( Int i = 0; i < sz; i++ ) cc += dwcopy[i][0] * dwcopy[i][0];

        cc -= l * l;
        Doub delta = bb * bb - 4. * aa * cc;
        Doub dlamb2 = ( -bb + sqrt ( delta ) ) / ( 2. * aa ); //maior
        Doub dlamb1= ( -bb - sqrt ( delta ) ) / ( 2. * aa ); //menor
        NRvector<Doub> v ( 2 );
        v[0]=dlamb1;
        v[1]=dlamb2;

        MatDoub temp1,sol1,temp2,sol2;
        temp1=dwb;
        temp1*=dlamb1;
        temp1+=dws;
        temp1+=dw;
        temp1.Mult ( dw,sol1 );

        return v;
    }



};

struct quadraticfunction   {
    Doub fa,fb,fc;
    Doub operator() ( const Doub x )
    {
        return fa*x*x+fb*x+fc;
    }
    Doub df ( const Doub x )
    {
        return 2*fa*x+fb;
    }
    void Set ( Doub a,Doub b,Doub c )
    {
        fa=a;
        fb=b;
        fc=c;
    }
};

