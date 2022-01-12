//
// Created by Diogo Cecílio on 10/12/21.
//

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


/**
 * @ingroup tools
 * @brief This class implements the necessary tools to solve a 3D elastic beam.
 */

using namespace Eigen;
class beam3dtools
{
public:
    /**
     * @brief Default constructor
     */
    beam3dtools();

    /**
     * @brief Class constructor
     * @param [in] mesh
     */
    beam3dtools ( mesh * inmesh );
    ~beam3dtools();



    /**
     * @brief Read a mesh from files
     * @param[in] filenameel file containing all the elements topology.
     * @param[in] filenamecoord file containing all the nodes coordinates.
     * @param[out] allcoords object containg the element topology and nodes in the format:
     * {{ {x1,y1,z1},{x2,y2,z2}...,{xn,yn,zn} }(el 1 topology), { {x1,y1,z1},{x2,y2,z2}...,{xn,yn,zn} }(el 2 topology),...,{ {x1,y1,z1},{x2,y2,z2}...{xn,yn,zn} (el m topology)}} n=1 to n=nelnode
     * @param[out] meshcoords object containg the mesh coordinates in the format:
     * { {x1,y1,z1},{x2,y2,z2}...,{xn,yn,zn} }  with  n=1 to n=nmeshnodes
     * @param[out] meshtopology object containg the mesh topology in the format:
     * { {id1,id2,id3...,idn}(el 1),{id1,id2,id3...,idn}(el 2)...,{id1,id2,id3...,idn}(el m) }  with  n=1 to n=nelnode
     */
    void ReadMesh ( std::vector<std::vector< std::vector<Doub > > >& allcoords, MatDoub& meshcoords, MatInt& meshtopology, string filenameel, string filenamecoord );

    /**
     * @brief Convert a string vector to a generic type T vector
     * @param[in] std::vector<std::string>  string vector.
     * @param[out] std::vector<T>  generic type T vector
     */
    template <class T>
    std::vector<T> vecstr_to_vec ( std::vector<std::string> vs );

    /**
     * @brief Get a specific element coordinates
     * @param[in] allcoords object containg the element topology and nodes.
     * @param[in] el  element that you want to get the coordinates
     * @param[out] elcoords  element coordinates
     */
    void GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, NRmatrix<Doub>  & elcoords );
    /**
     * @brief Solve a sparse linear system with Eigen
     * @param[in] A matrix A.
     * @param[in] b  vector b
     * @param[out] x  solution vector
     */
    void  SolveEigenSparse ( MatDoub A, MatDoub b, MatDoub& x );
    /**
     * @brief Solve a linear system with Eigen
     * @param[in] A matrix A.
     * @param[in] b  vector b
     * @param[out] x  solution vector
     */
    void SolveEigen ( MatDoub A, MatDoub b, MatDoub& x );
    /**
     * @brief Creates the material and the mesh for a 3D beam
     * @param[out] getmesh  finite element mesh
     * @param[out] mat  material
     */
    void CreateMatAndMesh ( mesh&getmesh, material &mat );
    /**
     * @brief Creates the material and the mesh for a cube
     * @param[out] getmesh  finite element mesh
     * @param[out] mat  material
     */
    void CreateMatAndMeshCube ( mesh&getmesh, material &mat );
    /**
     * @brief Manages other methods of the class and solves a elastic beam
     */
    void SolveElasticBeam();
    /**
     * @brief Manages other methods of the class and solves a elastic cube
     */
    void SolveElasticCube();

    /**
     * @brief Manages other methods of the class and solves a elastoplastic beam. Perform the newtom method for a nonlinear material
     */
    void IterativeProcess();

    /**
     * @brief Insert Dirichlet boundary conditions in the end of the beam (clamped in the end)
     * @param[in] K stiffness matrix.
     * @param[in] F  load vector
     * @param[in] mesh0  FEM mesh
     * @param[out] K modified stiffness matrix with BCs.
     * @param[out] F  modified load vector with BCs
     */
    void InsertBC ( mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F );

    /**
     * @brief Insert boundary conditions for the cube. Clamped in the base and a uniform pressure in the top
     * @param[in] K stiffness matrix.
     * @param[in] F  load vector
     * @param[in] mesh0  FEM mesh
     * @param[out] K modified stiffness matrix with BCs.
     * @param[out] F  modified load vector with BCs
     */
    void  InsertCubeBC ( mesh &mesh0,NRmatrix<Doub> & K, NRmatrix<Doub> & F );
    /**
     * @brief Compute a arc-length multiplier
     */
    Doub computelamda ( MatDoub& dwb, MatDoub& dws, MatDoub& dw, Doub& l );


    /**
     * @brief Search for a point, line or face ids.
     * @param[in] constcoord Constcoorddata vector containig info about the face you want to search for the id. If its a face, it must contain any coodinate locate in the in the face. If its a line, two coordinates must match. If its a point, tree exact coordinates must match
     * @param[in] allcoords  object containg the element topology and nodes
     * @param[in] meshtopology  object containg the mesh topology
     * @param[out] ids object(point, line or face) id or ids
     */
    void  FindIds ( NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<std::vector< std::vector<Doub > > >& allcoords, MatInt& meshtopology, std::vector<int>& ids );

    /**
     * @brief Insert boundary load BCs for the beam
     */
    void LoadBC ( mesh &mesh0, NRmatrix<Doub> & K, NRmatrix<Doub> & F );

    /**
     * @brief Write a file
     */
    void OutPutFile ( MatDoub & postdata, std::ofstream &file )
    {
        file.clear();
        for ( Int i = 0; i < postdata.nrows(); i++ ) {
            file << postdata[i][0] << " " << postdata[i][1] << endl;
        }

        file.close();
    }
    /**
     * @brief Write a file
     */
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


private:

    mesh* fmesh;

};
