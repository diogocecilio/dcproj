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
//#include "vtkmesh.h"
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

    void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords);
    void  SolveEigenSparse(MatDoub A, MatDoub b, MatDoub& x);
    void SolveEigen(MatDoub A, MatDoub b, MatDoub& x);
    void CreateMatAndMesh(mesh&getmesh, material &mat);
    void CreateMatAndMeshCube(mesh&getmesh, material &mat);
    void SolveElasticCube();
    void IterativeProcess();

    void InsertBC(mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F);

    void  InsertCubeBC(mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F);
    Doub computelamda(MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l);
    void  FindIds(NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& ids);
    void LoadBC(mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F);


        void OutPutFile(MatDoub & postdata, std::ofstream &file)
    {
        file.clear();
        for (Int i = 0;i < postdata.nrows(); i++)
        {
            file << postdata[i][0] << " " << postdata[i][1] << endl;
        }

        file.close();
    }

    void OutPutPost(std::vector<std::vector<double>> & postdata, std::ofstream &file)
    {
        file.clear();
        for (Int i = 0;i < postdata.size(); i++)
        {
            for (Int j = 0;j < postdata[0].size();j++)
            {
                file << postdata[i][j] << " ";
            }
            file << endl;
        }
        file.close();
    }


private:

    mesh* fmesh;

};
