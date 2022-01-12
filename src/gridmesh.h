#pragma once
#include "nr3.h"

class gridmesh
{
public:
    gridmesh ( Doub L, Doub h, Int nx, Int ny, Int order );
    ~gridmesh();

    void CreateMesh ( std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub & meshcoords, MatInt & meshtopology );
    void PrintAllCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, MatInt  meshtopology );
    static void GetElCoords ( std::vector<std::vector< std::vector<Doub > > > allcoords, Int el, MatDoub & elcoords );
    static void FindIdsInPath ( const MatDoub & path, std::vector<std::vector< std::vector<Doub > > > &allcoords, MatInt & meshtopology, std::vector<int> & idpath );
    static void Line ( VecDoub a, VecDoub b, Int ndivs, MatDoub & path );
    void PrintGMeshVTK ( std::vector<std::vector< std::vector<Doub > > >  allcoords, MatInt meshtopology, std::ofstream &file );

    inline void GetData ( std::vector<std::vector< std::vector<Doub > > > &allcoords, MatDoub  &meshcoords, MatInt &meshtopology )
    {
        allcoords = fallcoords;
        meshcoords = fmeshcoords;
        meshtopology = fmeshtopology;
    }


private:




public:

    std::vector<std::vector< std::vector<Doub > > > fallcoords;
    MatDoub  fmeshcoords;
    MatInt  fmeshtopology;
    Int forder;
    Doub fL;
    Doub fh;
    Int fnx;
    Int fny;

};

