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


class readgidmesh
{
public:
    /**
     * @brief Default constructor
     */
    readgidmesh();

    /**
     * @brief Class constructor
     * @param [in] mesh
     */
    readgidmesh ( string file );
    ~readgidmesh();




    template <class T>
    std::vector<T> str_vec( std::vector<std::string> &vs );


    void ReadMesh ( );
	

    
    void  FindIds ( NRvector<double> constcoorddata,NRvector<int> constcoord, std::vector<int>& ids )
    {
		//constcoorddata vector containig info about face to search id. It must contain any coodinate locate in the in the face
        MatDoub elcoords;
        int nels = fallcoords.size();
        GetElCoords (  0, elcoords );
        Int nnodes = elcoords.nrows();
        int sum=0;
        //constcoord.size() = 1 face
        //constcoord.size() = 2 linha
        //constcoord.size() = 3 pontos

        std::vector<int> dirs;
        for ( int iconst=0; iconst<constcoord.size(); iconst++ ) {
            sum+=constcoord[iconst];
            if ( constcoord[iconst]==1 ) {
                dirs.push_back ( iconst );
            }

        }
        for ( Int iel = 0; iel < nels; iel++ ) {
            GetElCoords (  iel, elcoords );
            for ( Int inode = 0; inode < nnodes; inode++ ) {

                if ( sum==1 ) {
                    if ( fabs ( elcoords[inode][dirs[0]] - constcoorddata[dirs[0]] ) <1.e-4 ) {
                        ids.push_back ( fmeshtopology[iel][inode] );
                    }
                } else if ( sum==2 ) {
                    if ( fabs ( elcoords[inode][dirs[0]] - constcoorddata[dirs[0]] ) <1.e-4 && abs ( elcoords[inode][dirs[1]] - constcoorddata[dirs[1]] ) <1.e-4 ) {
                        ids.push_back ( fmeshtopology[iel][inode] );
                    }
                } else if ( sum==3 ) {
                    if ( fabs ( elcoords[inode][dirs[0]] - constcoorddata[dirs[0]] ) <1.e-4 && abs ( elcoords[inode][dirs[1]] - constcoorddata[dirs[1]] ) <1.e-4 && abs ( elcoords[inode][dirs[2]] - constcoorddata[dirs[2]] ) <1.e-4 ) {
                        ids.push_back ( fmeshtopology[iel][inode] );
                    }
                }


            }


        }

        sort ( ids.begin(), ids.end() );
        ids.erase ( unique ( ids.begin(), ids.end() ), ids.end() );
    }
    void  GetElCoords ( Int el, NRmatrix<Doub>  & elcoords )
    {
        elcoords.assign ( fallcoords[el].size(), 3, 0. );
        for ( Int j = 0; j < fallcoords[el].size(); j++ ) {
            Doub x = fallcoords[el][j][0];
            Doub y = fallcoords[el][j][1];
            Doub z = fallcoords[el][j][2];
            elcoords[j][0] = x;
            elcoords[j][1] = y;
            elcoords[j][2] = z;
        }
    }
	
	inline NRmatrix<Int> GetTopology()
	{
		return fmeshtopology;
	}
	inline NRmatrix<Doub> GetCoords()
	{
		return fmeshcoords;
	}
	std::vector<std::vector< std::vector<Doub > > > GetAllCoords()
	{
		return fallcoords;
	}

    std::vector<std::vector<int>>   LineTopology ( std::vector<int> ids )
    {
        Int k = 0;


        std::vector<std::vector<int>> vg;
        for ( int j = 0; j < ids.size() / (fOrder+1); j++ ) {
            std::vector<int> v;
            for ( int i = 0; i < fOrder + 1; i++ ) {
                v.push_back ( ids[i + k] );
            }
            vg.push_back ( v );
            k += fOrder;
        }
        return vg;
    }

void  ToMatInt ( std::vector<std::vector<int>> in, MatInt & out )
{
    Int  rows = in.size();
    Int cols = in[0].size();
    out.assign ( rows, cols, 0. );
    for ( Int i = 0; i < rows; i++ ) for ( Int j = 0; j < cols; j++ ) out[i][j] = in[i][j];
}
	

private:

string ffile;
NRmatrix<Int> fmeshtopology;
NRmatrix<Doub> fmeshcoords;
std::vector<std::vector< std::vector<Doub > > > fallcoords;
Int fOrder;
};


