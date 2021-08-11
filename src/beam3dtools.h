#pragma once
#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"
#include <Eigen/Core>
#include <iostream>
#include "mesh.h"
#include "elastoplastic3D.h"
#include <chrono>
//using Eigen::MatrixXd;
using namespace Eigen;
class beam3dtools
{
public:

	beam3dtools();
    beam3dtools(mesh * inmesh);
	~beam3dtools();

    void SolveElasticBeam();
    void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord);
    template <class T>
    std::vector<T> vecstr_to_vec(std::vector<std::string> vs);
    void FindIdsInFace(NRvector<double> constcoorddata,int constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idsface);

    void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords);
    void  SolveEigenSparse(MatDoub A, MatDoub b, MatDoub& x);
    void SolveEigen(MatDoub A, MatDoub b, MatDoub& x);
    void CreateMatAndMesh(mesh&getmesh, material &mat);

    std::vector<std::vector<double>>   IterativeProcess(int ndesi, Doub dlamb0, Doub alphatol, int niter);

    void InsertBC(mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F);
    void InsertBC2(mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F);
    Doub computelamda(MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l);
private:

    mesh* fmesh;

};
