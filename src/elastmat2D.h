#pragma once

#include "shapequad.h"
#include "material.h"
#include "mesh.h"
#include "shapequad.h"
class elastmat2D
{
public:
    elastmat2D ( mesh &inmesh,Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order );
    elastmat2D ( Doub young, Doub nu, Doub thickness, Doub bodyforce, Int planestress, Int order, MatDoub  HHAT );
    elastmat2D();
    ~elastmat2D();

    void Contribute ( MatDoub &ek, MatDoub &ef, Doub xi, Doub eta, Doub w, MatDoub elcoords );
    void CalcStiff ( MatDoub &ek, MatDoub &ef, const MatDoub  &elcoords );
    //void Assemble(MatDoub &KG, MatDoub &FG, const std::vector<std::vector< std::vector<Doub > > > &allcoords, const MatDoub &meshnodes, const MatInt meshtopology);
    void assembleBandN ( MatDoub &B, MatDoub &N, const MatDoub &psis, const MatDoub &GradPhi );
    void assembleConstitutiveMatrix ( MatDoub &C, Doub mult );
    void GetElCoords ( std::vector<std::vector< std::vector<Doub > > > &allcoords, Int el, MatDoub & elcoords );
    void DirichletBC ( MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val );
    void ContributeLineNewan ( MatDoub &KG, MatDoub & FG, std::vector<int> ids, Int  dir, Int val );
    void SolPt ( const Int &el, const  MatDoub &solG, const Doub &xi, const Doub &eta, MatDoub &xycoords, MatDoub &sol );

    void PostProcess ( const MatDoub & nodalsol, std::vector<std::vector<double>> &solx, std::vector<std::vector<double>> &soly );
    void PostProcess ( const MatDoub & nodalsol, std::vector<std::vector<double>> &sol );


private:
    Doub fyoung;
    Doub fnu;
    Doub fbodyforce;
    Int fplanestress;
    Doub fthickness;
    Int fOrder;
    mesh fmesh;
public:
    MatDoub fHHAT;
    MatDoub fhhatvel;
};

