#pragma once
#pragma once


#include "nr3.h"



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
	slopeproject(mesh* inmesh,	KLGalerkinRF * inklgalerking);
	slopeproject(mesh* inmesh, KLGalerkinRF* inklgalerking, NRmatrix<MatDoub> randomfield);
	slopeproject(mesh inmesh, KLGalerkinRF inklgalerking, NRmatrix<MatDoub> randomfield);
	slopeproject();
	~slopeproject();
	slopeproject(slopeproject& copy);
	void CreateRandomField(string namefolder);
    void PrintMCS(string namefolder,int imc,bool print);
	void OutPutPost(MatDoub& postdata, std::ofstream& file);
	void OutPutPost(MatInt& postdata, std::ofstream& file);
	void OutPutPost(std::vector<std::vector<double>>& postdata, std::ofstream& file);
	void OutPutFile(MatDoub& postdata, std::ofstream& file);
	void OutPutFile1var(MatDoub& postdata, std::ofstream& file);
	void OutPutFile4var(MatDoub& postdata, std::ofstream& file);
    void PrintMathematicaFormat(MatDoub postdata, std::ofstream& file);
	Doub computelamda(MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l);
    void OutPutPost2(std::vector<std::vector<double>>& postdata, std::ofstream& file);
	void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord);
	std::vector<Doub>   vecstr_to_vecdoub2(std::vector<string> vs);
	template <class T>
	std::vector<T>   vecstr_to_vec(std::vector<string> vs);
	void  ReadMatDoub(MatDoub& matdoub, std::string  file);
	std::vector<Doub>   vecstr_to_vecdoub(std::vector<string> vs);
	std::vector<Int>   vecstr_to_vecint(std::vector<string> vs);
	std::vector<std::vector<double>>   IterativeProcess(int ndesi, Doub dlamb0, Doub alphatol, int niter);
	std::vector<std::vector<double>>   IterativeProcessArcLengthSRM(int ndesi, Doub dlamb0, Doub alphatol, int niter);
	void   IterativeProcess2();
    std::vector<std::vector<double>>    IterativeProcessGIMBinarySearch();
    void   IterativeProcessSRMBinarySearch();
    void SolveEigenSparse(MatDoub A, MatDoub b, MatDoub& x);
    void SolveEigen(SparseMatrix<double> A, VectorXd b, VectorXd& x);
	void SolveEigen(MatDoub A, MatDoub b, MatDoub& x);
    void SolveEigen2(MatDoub A, MatDoub b, MatDoub& x);
    void SolveEigen3(MatDoub A, MatDoub b, MatDoub& x);
	void InserBC(MatDoub& KG, MatDoub& R, MatDoub& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material* mat);
    void InserBC(SparseMatrix<double> & KG, VectorXd& R, VectorXd& FBODY, std::vector<int> idsbottom, std::vector<int> idsright, std::vector<int> idsleft, material*mat);
	void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub& elcoords);
	void FindIdsInPath(const MatDoub& path, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idpath);
	void Line(VecDoub a, VecDoub b, Int ndivs, MatDoub& path);

	void findbcids(mesh* gmesh, std::vector<std::vector<int>>& idsvector);
	MatDoub  AssembleHhationho(Int i);
	std::vector<std::vector<double>>  IterativeProcessShearRed(Doub fac, Doub delta,Doub tol);
    std::vector<std::vector<double>>  IterativeProcessShearRed2(Doub fac, Doub delta,Doub tol);
	void MonteCarloSRM(int iter, int iter2, bool print, string writenamefolder);
	void MonteCarloGIM(int iter, int iter2,bool print, string writenamefolder);
	mesh* fmesh;
	KLGalerkinRF* fklgalerking;
	NRmatrix<MatDoub> frandomfield;

};
