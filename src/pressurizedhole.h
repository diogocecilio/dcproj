#pragma once

#include "nr3.h"
#include <math.h>
#include <cmath>
#include "elastoplasticbase.h"
#include <Eigen/Core>
#include <iostream>
#include "mesh.h"
#include "elastoplastic2D.h"
#include <chrono>
#include "gridmesh.h"

//using Eigen::MatrixXd;
using namespace Eigen;
class pressurizedhole
{
public:

	pressurizedhole();
    pressurizedhole(mesh * inmesh);
	~pressurizedhole();

    void SolveElasticHole();
    void ReadMesh(std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord);
    template <class T>
    std::vector<T> vecstr_to_vec(std::vector<std::string> vs);
    void FindIdsInFace(NRvector<double> constcoorddata,int constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& idsface);

    void GetElCoords(std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords);
    void  SolveEigenSparse(MatDoub A, MatDoub b, MatDoub& x);
    void SolveEigen(MatDoub A, MatDoub b, MatDoub& x);
    void CreateMatAndMesh(mesh&getmesh, material &mat);

    void   IterativeProcess();

    void InsertBC(mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F);
    void LoadBC(mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F);
    Doub computelamda(MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l);
    std::vector<std::vector<int>>  LineTopology(std::vector<int> ids, Int order);
    void ToMatInt(std::vector<std::vector<int>> in, MatInt & out);
    //constcoorddata specify a coordinate on the element
    //constcoord[3] // specify with 1 in the fixed direction
     void  FindIds(NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& ids);


    template <class T>
    void OutPutPost(NRmatrix<T>& postdata, std::ofstream& file)
    {
        file.clear();
        for (Int i = 0; i < postdata.nrows(); i++)
        {
            for (Int j = 0; j < postdata.ncols(); j++)
            {
                file << postdata[i][j] << " ";
            }
            file << endl;
        }
        file.close();
    }

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
